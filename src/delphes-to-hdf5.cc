#include "delphes_utils.cc"
#include "gen_utils.cc"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <bitset>
#include <cstdint>
#include <highfive/highfive.hpp>

#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <indicators/progress_bar.hpp>

#include <filesystem>

namespace fs = std::filesystem;
using namespace HighFive;
using namespace indicators;


static const int16_t kElectronType = 0;
static const int16_t kMuonType = 1;


struct METData {
    float MET;
    float phi;
};

struct JetData {
    float pt;
    float eta;
    float phi;
    float energy;
    int8_t is_tagged;
};

struct LeptonData {
    float pt;
    float eta;
    float phi;
    float energy;
    int16_t charge;
    int16_t type;
};

struct ParticleData {
    int32_t PDGID;
    float pt;
    float eta;
    float phi;
    float mass;
};

inline CompoundType create_compound_MET() {
    return {{"MET", create_datatype<float>()},
            {"phi", create_datatype<float>()}};
}

inline CompoundType create_compound_Jet() {
    return {{"pt", create_datatype<float>()},
            {"eta", create_datatype<float>()},
            {"phi", create_datatype<float>()},
            {"energy", create_datatype<float>()},
            {"is_tagged", create_datatype<int8_t>()}};
}

inline CompoundType create_compound_Lepton() {
    return {{"pt", create_datatype<float>()},
            {"eta", create_datatype<float>()},
            {"phi", create_datatype<float>()},
            {"energy", create_datatype<float>()},
            {"charge", create_datatype<int16_t>()},
            {"type", create_datatype<int16_t>()}};
}

inline CompoundType create_compound_Particle() {
    return {{"PDGID", create_datatype<int32_t>()},
            {"pt", create_datatype<float>()},
            {"eta", create_datatype<float>()},
            {"phi", create_datatype<float>()},
            {"mass", create_datatype<float>()}};
}

// Register compound types with HighFive
HIGHFIVE_REGISTER_TYPE(METData, create_compound_MET)
HIGHFIVE_REGISTER_TYPE(JetData, create_compound_Jet)
HIGHFIVE_REGISTER_TYPE(LeptonData, create_compound_Lepton)
HIGHFIVE_REGISTER_TYPE(ParticleData, create_compound_Particle)


// Returns top and anti-top quarks
// Both particles are last copies
// Throws if multiple top or anti-top quarks are found or if one of them is missing
std::pair<const GenParticle*, const GenParticle*>
findTT(const std::vector<const GenParticle*>& particle_vec) {
    const GenParticle* t0 = nullptr; // top quark
    const GenParticle* t1 = nullptr; // anti-top quark
  for (auto particle : particle_vec) {
    if (particle->PID == 6) {
      // FIXME: hardcoded top quark PID
      particle = getLastCopy(particle, particle_vec);

      if (t0 == nullptr) {
        t0 = particle;

      } else if (t0 != particle) {
        throw std::runtime_error("found multiple top quarks in the event");

      }

    } else if (particle->PID == -6) {
      // FIXME: hardcoded anti top quark PID
      particle = getLastCopy(particle, particle_vec);

      if (t1 == nullptr) {
        t1 = particle;

      } else if (t1 != particle) {
        throw std::runtime_error("found multiple anti-top quarks in the event");

      }
    }

    if ((t0 != nullptr) and (t1 != nullptr)) {
      break;
    }
  }

  if ((t0 == nullptr) or (t1 == nullptr)) {
    throw std::runtime_error("failed to find top and anti-top quarks in the event");
  }
  return { t0, t1 };
}


// Returns two daughters of a mother particle
// A mother particle might have more than two daughters when the decay is
// simulated with Pythia8
std::pair<const GenParticle*, const GenParticle*>
findTwoDaughters(const GenParticle* mother,
                 const std::vector<const GenParticle*>& p_vec)
{

  const auto daughters = getDaughters(mother, p_vec);
  if (daughters.size() != 2) {
    throw std::runtime_error("expected exactly two daughters for particle PID " +
                             std::to_string(mother->PID) + ", found " +
                             std::to_string(daughters.size()));
  }

  const GenParticle* d1 = daughters.at(0);
  const GenParticle* d2 = daughters.at(1);

  return { d1, d2 };
}

// Returns two daughters of top: (b, W)
std::pair<const GenParticle*, const GenParticle*>
findTopQuarkDaughters(const GenParticle* mother,
                      const std::vector<const GenParticle*>& p_vec)
{
  auto [b, w] = findTwoDaughters(mother, p_vec);
  // swap if necessary
  if (std::abs(w->PID) != 24) { // FIXME: hardcoded W PID
    std::swap(b, w);
  }
  // b = getFirstCopy(b, p_vec);
  w = getLastCopy(w, p_vec);
  return { b, w };
}


// Returns two daughters of W: (lepton, neutrino) or (down-type quark, up-type quark)
// Both daughters are first copies
std::pair<const GenParticle*, const GenParticle*>
findWBosonDaughters(const GenParticle* mother,
                    const std::vector<const GenParticle*>& p_vec)
{
  auto [d1, d2] = findTwoDaughters(mother, p_vec);
  // swap if necessary to have lepton first if applicable
  // hadronic decay: down-type quark, up-type quark
  // leptonic decay: lepton, neutrino
  if (std::abs(d1->PID) > std::abs(d2->PID)) {
    std::swap(d1, d2);
  }

  return { d1, d2 };
}


// Returns a boolean indicating whether the specified bit is set in the value
bool
checkBit(const int value, const int bit)
{
  return (value & (1 << bit)) != 0;
}


// Main processing function
void run(
  const fs::path& input_path,
  const fs::path& output_path,
  const long max_events,
  const bool enable_progress_bar) {

  spdlog::info("Setting up Delphes...");
  setupDelphes();
  spdlog::info("Delphes setup complete.");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE: Input Delphes ROOT file and tree
  //
  //////////////////////////////////////////////////////////////////////////////

  spdlog::info("Opening input file: {}", input_path.string());
  auto input_file = TFile::Open(input_path.c_str());
  auto input_tree = dynamic_cast<TTree*>(input_file->Get("Delphes"));

  TClonesArray* input_jet_array = nullptr;
  TClonesArray* input_electron_array = nullptr;
  TClonesArray* input_muon_array = nullptr;
  TClonesArray* input_met_array = nullptr;
  TClonesArray* input_particle_array = nullptr;

  input_tree->SetBranchAddress("Jet", &input_jet_array);
  input_tree->SetBranchAddress("Electron", &input_electron_array);
  input_tree->SetBranchAddress("Muon", &input_muon_array);
  input_tree->SetBranchAddress("MissingET", &input_met_array);
  input_tree->SetBranchAddress("Particle", &input_particle_array);

  spdlog::info("Input file and tree setup complete.");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE: Output HDF5 file and datasets
  //
  //////////////////////////////////////////////////////////////////////////////

  spdlog::info("Creating output file: {}", output_path.string());

  const size_t output_chunk_size = 1024;
  const size_t max_jets_size = 10;
  const size_t lepton_size = 2; // l- and l+ (particle and anti-particle)
  const size_t neutrino_size = 2; // nu and nu~
  const size_t truth_particles_size = 10; // t, b, W, l, nu for both top and anti-top
  const size_t truth_quarks_size = 2; // b and b~

  spdlog::info("  - output chunk size: {}", output_chunk_size);
  spdlog::info("  - max jets size: {}", max_jets_size);
  spdlog::info("  - lepton size: {}", lepton_size);
  spdlog::info("  - neutrino size: {}", neutrino_size);
  spdlog::info("  - truth particles size: {}", truth_particles_size);
  spdlog::info("  - truth quarks size: {}", truth_quarks_size);

  HighFive::File output_file{output_path.string(), HighFive::File::Truncate};
  HighFive::Group output_group = output_file.createGroup("delphes");

  // met dataset
  spdlog::info("  - Creating MET dataset...");
  DataSpace met_space({0}, {DataSpace::UNLIMITED});
  DataSetCreateProps met_props{};
  met_props.add(Chunking{std::vector<hsize_t>{output_chunk_size}});
  auto met_dataset = output_group.createDataSet<METData>("MET", met_space, met_props);

  // decay_channel
  spdlog::info("  - Creating decay_channel dataset...");
  DataSpace decay_channel_space({0}, {DataSpace::UNLIMITED});
  DataSetCreateProps decay_channel_props{};
  decay_channel_props.add(Chunking{std::vector<hsize_t>{output_chunk_size}});
  auto decay_channel_dataset = output_group.createDataSet<int16_t>("decay_channel", decay_channel_space, decay_channel_props);

  // jets
  spdlog::info("  - Creating jets dataset...");
  DataSpace jets_space({0, max_jets_size}, {DataSpace::UNLIMITED, max_jets_size});
  DataSetCreateProps jets_props{};
  jets_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, max_jets_size}});
  auto jets_dataset = output_group.createDataSet<JetData>("jets", jets_space, jets_props);

  // jets_indices
  spdlog::info("  - Creating jets_indices dataset...");
  DataSpace jets_indices_space({0, max_jets_size}, {DataSpace::UNLIMITED, max_jets_size});
  DataSetCreateProps jets_indices_props{};
  jets_indices_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, max_jets_size}});
  auto jets_indices_dataset = output_group.createDataSet<int16_t>("jets_indices", jets_indices_space, jets_indices_props);

  // leptons
  spdlog::info("  - Creating leptons dataset...");
  DataSpace leptons_space({0, lepton_size}, {DataSpace::UNLIMITED, lepton_size});
  DataSetCreateProps leptons_props;
  leptons_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, lepton_size}});
  auto leptons_dataset = output_group.createDataSet<LeptonData>("leptons", leptons_space, leptons_props);

  // matchability
  spdlog::info("  - Creating matchability dataset...");
  DataSpace matchability_space({0}, {DataSpace::UNLIMITED});
  DataSetCreateProps matchability_props;
  matchability_props.add(Chunking{std::vector<hsize_t>{output_chunk_size}});
  auto matchability_dataset = output_group.createDataSet<int16_t>("matchability", matchability_space, matchability_props);

  // nbjets
  spdlog::info("  - Creating nbjets dataset...");
  DataSpace nbjets_space({0}, {DataSpace::UNLIMITED});
  DataSetCreateProps nbjets_props;
  nbjets_props.add(Chunking{std::vector<hsize_t>{output_chunk_size}});
  auto nbjets_dataset = output_group.createDataSet<int16_t>("nbjets", nbjets_space, nbjets_props);

  // neutrinos
  spdlog::info("  - Creating neutrinos dataset...");
  DataSpace neutrinos_space({0, neutrino_size}, {DataSpace::UNLIMITED, neutrino_size});
  DataSetCreateProps neutrinos_props{};
  neutrinos_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, neutrino_size}});
  auto neutrinos_dataset = output_group.createDataSet<ParticleData>("neutrinos", neutrinos_space, neutrinos_props);

  // njets
  spdlog::info("  - Creating njets dataset...");
  DataSpace njets_space({0}, {DataSpace::UNLIMITED});
  DataSetCreateProps njets_props;
  njets_props.add(Chunking{std::vector<hsize_t>{output_chunk_size}});
  auto njets_dataset = output_group.createDataSet<int16_t>("njets", njets_space, njets_props);

  // truth_particles
  spdlog::info("  - Creating truth_particles dataset...");
  DataSpace truth_particles_space({0, truth_particles_size}, {DataSpace::UNLIMITED, truth_particles_size});
  DataSetCreateProps truth_particles_props{};
  truth_particles_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, truth_particles_size}});
  auto truth_particles_dataset = output_group.createDataSet<ParticleData>("truth_particles", truth_particles_space, truth_particles_props);

  // truth_quarks
  spdlog::info("  - Creating truth_quarks dataset...");
  DataSpace truth_quarks_space({0, truth_quarks_size}, {DataSpace::UNLIMITED, truth_quarks_size});
  DataSetCreateProps truth_quarks_props{};
  truth_quarks_props.add(Chunking{std::vector<hsize_t>{output_chunk_size, truth_quarks_size}});
  auto truth_quarks_dataset = output_group.createDataSet<ParticleData>("truth_quarks", truth_quarks_space, truth_quarks_props);

  spdlog::info("Output file and datasets setup complete.");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE: Event selection criteria
  //
  // TODO: move to config file?
  //
  //////////////////////////////////////////////////////////////////////////////

  const float jet_pt_min = 25.0;
  const float jet_abs_eta_max = 2.5;
  const uint64_t jet_btag_bit = 0; // NOTE: https://github.com/delphes/delphes/blob/master/cards/delphes_card_CMS.tcl#L720
  const float electron_pt_min = 15.0;
  const float electron_abs_eta_max = 2.5;
  const float muon_pt_min = 15.0;
  const float muon_abs_eta_max = 2.5;
  const float parton_jet_matching_max_dr = 0.4;

  spdlog::info("Event selection criteria:");
  spdlog::info("  - jet pT > {} GeV and |eta| < {}", jet_pt_min, jet_abs_eta_max);
  spdlog::info("  - electron pT > {} GeV and |eta| < {}", electron_pt_min, electron_abs_eta_max);
  spdlog::info("  - muon pT > {} GeV and |eta| < {}", muon_pt_min, muon_abs_eta_max);
  spdlog::info("  - parton-jet matching max deltaR = {}", parton_jet_matching_max_dr);


  //////////////////////////////////////////////////////////////////////////////
  // NOTE: Event loop
  //
  //////////////////////////////////////////////////////////////////////////////

  ProgressBar progress_bar{};
  if (enable_progress_bar) {
    spdlog::debug("Enabling progress bar...");

    progress_bar.set_option(option::BarWidth{ 80 });
    progress_bar.set_option(option::Start{ " [" });
    progress_bar.set_option(option::Fill{ "=" });
    progress_bar.set_option(option::Lead{ ">" });
    progress_bar.set_option(option::Remainder{ "-" });
    progress_bar.set_option(option::End{ "]" });
    // progress_bar.set_option(option::PrefixText{ "Processing... 👀" });
    progress_bar.set_option(option::ForegroundColor{ Color::green });
    progress_bar.set_option(
      option::FontStyles{ std::vector<FontStyle>{ FontStyle::bold } });
    progress_bar.set_option(option::ShowPercentage{ true });
    progress_bar.set_option(option::ShowElapsedTime{ true });
    progress_bar.set_option(option::ShowRemainingTime{ true });

  } else {
    spdlog::debug("Progress bar disabled.");

  }

  const auto total = static_cast<size_t>(input_tree->GetEntries());
  const auto entry_stop = (max_events > 0) ? std::min(static_cast<size_t>(max_events), total) : total;

  size_t num_selected = 0;
  // TODO: cutflow?

  spdlog::info("Starting event loop...");
  spdlog::info("  - Total entries in input file: {}", total);
  spdlog::info("  - Entries to process: {}", entry_stop);

  for (size_t entry = 0; entry < entry_stop; ++entry) {
    input_tree->GetEntry(entry);

    if (enable_progress_bar) {
      progress_bar.set_progress(100 * (static_cast<double>(entry + 1) / total));
    }

    // NOTE: jet selection
    // At least two jets with pT > 25 GeV in the range |eta|<2.5 are required.
    std::vector<const Jet*> jet_vec{};

    const size_t jet_stop = std::min(static_cast<size_t>(input_jet_array->GetEntries()), max_jets_size);
    for (size_t jet_idx = 0; jet_idx < jet_stop; ++jet_idx) {
      if (const auto jet = dynamic_cast<const Jet*>(input_jet_array->At(jet_idx))) {
        if ((jet->PT > jet_pt_min) and (std::abs(jet->Eta) < jet_abs_eta_max)) {
          jet_vec.push_back(jet);
        }

      } else {
        throw std::runtime_error("failed to dynamic_cast Jet");

      }
    }
    jet_vec.shrink_to_fit();

    // NOTE: NOTE: lepton selection
    std::vector<std::variant<const Electron*, const Muon*> > lepton_vec{};

    // electrons
    for (int idx = 0; idx < input_electron_array->GetEntries(); ++idx) {
      if (const auto electron = dynamic_cast<const Electron*>(input_electron_array->At(idx))) {
        if ((electron->PT > electron_pt_min) and (std::abs(electron->Eta) < electron_abs_eta_max)) {
          lepton_vec.push_back(electron);
        }
      } else {
        throw std::runtime_error("failed to dynamic_cast Electron");
      }

    }
    // muons
    for (int idx = 0; idx < input_muon_array->GetEntries(); ++idx) {
      if (const auto muon = dynamic_cast<const Muon*>(input_muon_array->At(idx))) {
        if ((muon->PT > muon_pt_min) and (std::abs(muon->Eta) < muon_abs_eta_max)) {
          lepton_vec.push_back(muon);
        }
      } else {
        throw std::runtime_error("failed to dynamic_cast Muon");
      }
    }
    lepton_vec.shrink_to_fit();

    // NOTE: missing transverse momentum
    const auto met = dynamic_cast<const MissingET*>(input_met_array->At(0));

    ////////////////////////////////////////////////////////////////////////////
    // NOTE: event selection
    //
    ////////////////////////////////////////////////////////////////////////////

    // NOTE: at least two jets
    if (jet_vec.size() < 2) {
      continue;
    }

    // NOTE: exactly two leptons
    if (lepton_vec.size() != 2) {
      continue;
    }

    // NOTE: leptons must have opposite charge
    auto get_lepton_charge = [](auto&& arg) -> int {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, const Electron*>) {
          return arg->Charge;

        } else if constexpr (std::is_same_v<T, const Muon*>) {
          return arg->Charge;

        } else {
          throw std::runtime_error("unknown variant type");

        }
    };

    const auto l0_charge = std::visit(get_lepton_charge, lepton_vec.at(0));
    const auto l1_charge = std::visit(get_lepton_charge, lepton_vec.at(1));

    if (l0_charge * l1_charge >= 0) {
      continue;
    }

    ////////////////////////////////////////////////////////////////////////////
    // NOTE: GenParticle
    //
    ////////////////////////////////////////////////////////////////////////////
    std::vector<const GenParticle*> gen_particle_vec{};
    for (int idx = 0; idx < input_particle_array->GetEntries(); ++idx) {
      const auto particle = dynamic_cast<const GenParticle*>(input_particle_array->At(idx));
      gen_particle_vec.push_back(particle);
    }

    // t0: t, t1: t~
    const auto [t0, t1] = findTT(gen_particle_vec);
    // b0: b, w0: W+
    auto [b0, w0] = findTopQuarkDaughters(t0, gen_particle_vec);
    // b1: b~, w1: W-
    auto [b1, w1] = findTopQuarkDaughters(t1, gen_particle_vec);

    w0 = getLastCopy(w0, gen_particle_vec);
    w1 = getLastCopy(w1, gen_particle_vec);

    // lep_0: l+, nu_0: vl
    auto [lep_0, nu_0] = findWBosonDaughters(w0, gen_particle_vec);
    // lep_1: l-, nu_1: vl~
    auto [lep_1, nu_1] = findWBosonDaughters(w1, gen_particle_vec);

    // NOTE: parton-jet matching

    // b and b~
    const std::vector<const GenParticle*> parton_vec{b0, b1};
    const auto parton_jet_idx_vec =
      matchPartonToJet(parton_vec, jet_vec, /*max_distance=*/parton_jet_matching_max_dr);

    ////////////////////////////////////////////////////////////////////////////
    // NOTE: create data to write
    //
    ////////////////////////////////////////////////////////////////////////////

    // met
    const METData met_data{ met->MET, met->Phi };

    // decay channel
    // FIXME: assume both tops decay leptonically
    // This should be inferred from the decay products of the W bosons
    std::bitset<2> decay_channel_bits{};
    decay_channel_bits.set(0, true); // t0 decays leptonically
    decay_channel_bits.set(1, true); // t1 decays leptonically

    const int16_t decay_channel_data = static_cast<int16_t>(decay_channel_bits.to_ulong());

    // jets data
    std::vector<JetData> jets_data{};
    jets_data.reserve(max_jets_size);
    for (const auto& jet : jet_vec) {
      jets_data.emplace_back(
        jet->PT,
        jet->Eta,
        jet->Phi,
        static_cast<float>(jet->P4().Energy()),
        static_cast<int8_t>(checkBit(jet->BTag, jet_btag_bit) ? 1 : 0)
      );
    }
    // pad with zeros
    for (size_t idx = jet_vec.size(); idx < max_jets_size; ++idx) {
      jets_data.emplace_back(0.0f, 0.0f, 0.0f, 0.0f, 0);
    }

    // jets_indices
    std::vector<int16_t> jets_indices_data(jet_vec.size(), -1);
    for (size_t parton_idx = 0; parton_idx < parton_vec.size(); ++parton_idx) {
      const auto& jet_idx = parton_jet_idx_vec.at(parton_idx);
      if (jet_idx == -1) {
        continue;
      }
      jets_indices_data.at(jet_idx) = static_cast<int16_t>(parton_idx);
    }
    // pad with -1
    for (size_t idx = jet_vec.size(); idx < max_jets_size; ++idx) {
      jets_indices_data.emplace_back(-1);
    }

    // leptons
    std::vector<LeptonData> leptons_data{};
    leptons_data.reserve(lepton_size);
    for (const auto& lepton : lepton_vec) {
      std::visit([&leptons_data](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, const Electron*>) {
          leptons_data.emplace_back(
            arg->PT,
            arg->Eta,
            arg->Phi,
            static_cast<float>(arg->P4().Energy()),
            static_cast<int16_t>(arg->Charge),
            kElectronType
          );

        } else if constexpr (std::is_same_v<T, const Muon*>) {
          leptons_data.emplace_back(
            arg->PT,
            arg->Eta,
            arg->Phi,
            static_cast<float>(arg->P4().Energy()),
            static_cast<int16_t>(arg->Charge),
            kMuonType
          );

        } else {
          throw std::runtime_error("unknown variant type");

        }
      }, lepton);
    }

    // NOTE: matchability
    std::bitset<2> matchability_bits{};
    for (size_t parton_idx = 0; parton_idx < parton_vec.size(); ++parton_idx) {
      matchability_bits.set(parton_idx, parton_jet_idx_vec.at(parton_idx) >= 0); // b
    }
    const int16_t matchability_data = static_cast<int16_t>(matchability_bits.to_ulong());

    // NOTE: nbjets
    const int16_t nbjets_data = static_cast<int16_t>(
      std::count_if(
        jet_vec.cbegin(),
        jet_vec.cend(),
        [](const Jet* jet) { return checkBit(jet->BTag, jet_btag_bit); }
      )
    );

    // NOTE: neutrinos data
    std::vector<ParticleData> neutrinos_data{};
    neutrinos_data.reserve(neutrino_size);

    for (const auto& nu : {nu_0, nu_1}) {
      neutrinos_data.emplace_back(
        static_cast<int32_t>(nu->PID),
        nu->PT,
        nu->Eta,
        nu->Phi,
        static_cast<float>(nu->Mass)
      );
    }

    // NOTE: njets
    const auto njets_data = static_cast<int16_t>(jet_vec.size());

    // NOTE:
    // truth particles data
    std::vector<ParticleData> truth_particles_data{};
    truth_particles_data.reserve(truth_particles_size);

    const auto truth_particle_vec = {
      t0, b0, w0, lep_0, nu_0,
      t1, b1, w1, lep_1, nu_1,
    };

    for (const auto& particle : truth_particle_vec) {
      truth_particles_data.emplace_back(
        static_cast<int32_t>(particle->PID),
        particle->PT,
        particle->Eta,
        particle->Phi,
        static_cast<float>(particle->Mass)
      );
    }

    // NOTE: truth quarks data
    std::vector<ParticleData> truth_quarks_data{};
    truth_quarks_data.reserve(truth_quarks_size);

    for (const auto& quark : {b0, b1}) {
      truth_quarks_data.emplace_back(
        static_cast<int32_t>(quark->PID),
        quark->PT,
        quark->Eta,
        quark->Phi,
        static_cast<float>(quark->Mass)
      );
    }

    ////////////////////////////////////////////////////////////////////////////
    // NOTE:: resize datasets
    //
    ////////////////////////////////////////////////////////////////////////////
    num_selected += 1;
    // index of the newly added entry
    const size_t output_idx = num_selected - 1;

    met_dataset.resize({num_selected});
    decay_channel_dataset.resize({num_selected});
    jets_dataset.resize({num_selected, max_jets_size});
    jets_indices_dataset.resize({num_selected, max_jets_size});
    leptons_dataset.resize({num_selected, lepton_size});
    matchability_dataset.resize({num_selected});
    nbjets_dataset.resize({num_selected});
    neutrinos_dataset.resize({num_selected, neutrino_size});
    njets_dataset.resize({num_selected});
    truth_particles_dataset.resize({num_selected, truth_particles_size});
    truth_quarks_dataset.resize({num_selected, truth_quarks_size});

    // NOTE: write data
    met_dataset.select({output_idx}, {1}).write(&met_data);
    decay_channel_dataset.select({output_idx}, {1}).write(&decay_channel_data);
    jets_dataset.select({output_idx, 0}, {1, max_jets_size}).write_raw(jets_data.data());
    jets_indices_dataset.select({output_idx, 0}, {1, max_jets_size}).write_raw(jets_indices_data.data());
    leptons_dataset.select({output_idx, 0}, {1, lepton_size}).write_raw(leptons_data.data());
    matchability_dataset.select({output_idx}, {1}).write(&matchability_data);
    nbjets_dataset.select({output_idx}, {1}).write(&nbjets_data);
    neutrinos_dataset.select({output_idx, 0}, {1, neutrino_size}).write_raw(neutrinos_data.data());
    njets_dataset.select({output_idx}, {1}).write(&njets_data);
    truth_particles_dataset.select({output_idx, 0}, {1, truth_particles_size}).write_raw(truth_particles_data.data());
    truth_quarks_dataset.select({output_idx, 0}, {1, truth_quarks_size}).write_raw(truth_quarks_data.data());

  }

  const auto unweighted_efficiency =
    static_cast<double>(num_selected) / static_cast<double>(entry_stop);

  spdlog::info("Event loop complete.");
  spdlog::info("  - Processed entries: {}", entry_stop);
  spdlog::info("  - Selected entries: {}", num_selected);
  spdlog::info("  - Unweighted selection efficiency: {:.2f}%", 100 * unweighted_efficiency);
}


int
main(int argc, char* argv[])
{
  spdlog::info("delphes-to-hdf5 started");

  //////////////////////////////////////////////////////////////////////////////
  /// NOTE: Argument parser
  ///
  //////////////////////////////////////////////////////////////////////////////
  argparse::ArgumentParser parser{ "delphes-to-hdf5" };

  parser.add_argument("-i", "--input")
    .required()
    .help("input delphes root file");

  parser.add_argument("-o", "--output").required().help("output hdf5 file");

  parser.add_argument("--max-events")
    .default_value(-1l)
    .scan<'i', long>()
    .help(
      "maximum number of events to analyze (default: -1, analyze all events)");

  parser.add_argument("--progress")
    .default_value(false)
    .implicit_value(true)
    .help("enable progress bar");

  try {
    parser.parse_args(argc, argv);

  } catch (const std::exception& err) {
    spdlog::error("{}", err.what());
    std::cerr << parser << std::endl;
    return 1;
  }

  const fs::path input_path = parser.get<std::string>("input");
  const fs::path output_path = parser.get<std::string>("output");
  const long max_events = parser.get<long>("max-events");
  const bool enable_progress_bar = parser.get<bool>("progress");

  spdlog::info("Arguments:");
  spdlog::info("  - input file: {}", input_path.string());
  spdlog::info("  - output file: {}", output_path.string());
  spdlog::info("  - max events: {}", max_events);
  spdlog::info("  - progress bar: {}", enable_progress_bar ? "enabled" : "disabled");

  try {
    if (not fs::exists(input_path)) {
      throw std::runtime_error("input file not found: " + input_path.string());
    }

    if (fs::exists(output_path)) {
      throw std::runtime_error("output file already exists: " +
                               output_path.string());
    }

    run(
      input_path, output_path, max_events, enable_progress_bar);

  } catch (const std::exception& err) {
    spdlog::error("{}", err.what());
    return 1;

  }

  spdlog::info("delphes-to-hdf5 finished successfully.");
  return 0;
}
