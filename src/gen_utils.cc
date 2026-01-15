#include "Math/GenVector/VectorUtil.h"
#include "classes/DelphesClasses.h"

// include spdlog for logging
#include <spdlog/spdlog.h>

#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

using namespace ROOT::Math;

std::vector<const GenParticle*>
getDaughters(const GenParticle* mother,
             const std::vector<const GenParticle*>& p_arr);

const GenParticle*
getHardProcessMotherCopy(const GenParticle*,
                         const std::vector<const GenParticle*>&);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L511-L523
const GenParticle*
getPreviousCopy(const GenParticle* p,
                const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L525-L537
const GenParticle*
getNextCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L392-L404
const GenParticle*
getFirstCopy(const GenParticle* p,
             const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L406-L418
const GenParticle*
getLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

bool
isFirstCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

bool
isLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L265-L295
// TODO
bool
isHardProcess(const GenParticle* p);

bool
isFromHardProcess(const GenParticle*, const std::vector<const GenParticle*>&);

const GenParticle*
findDaughter(const GenParticle* mother,
             const int pid,
             const std::vector<const GenParticle*>& p_arr);

const GenParticle*
findLastDaughter(const GenParticle* mother,
                 const int pid,
                 const std::vector<const GenParticle*>& p_arr);

// Returns jet indices matched to partons
// match parton to jet with shortest distance
// starting with the shortest distance available
// apply some outlier rejection if desired
// https://github.com/cms-sw/cmssw/blob/CMSSW_14_1_0_pre4/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L149-L196
std::vector<int32_t>
matchPartonToJet(const std::vector<const GenParticle*> parton_vec,
                 const std::vector<const Jet*> jet_vec,
                 const float max_distance);

std::pair<const GenParticle*, const GenParticle*>
findTwoDaughters(const GenParticle*, const std::vector<const GenParticle*>&);

// Returns b and W
std::pair<const GenParticle*, const GenParticle*>
findTopQuarkDaughters(const GenParticle*,
                      const std::vector<const GenParticle*>&);

std::pair<const GenParticle*, const GenParticle*>
findWBosonDaughters(const GenParticle*, const std::vector<const GenParticle*>&);

std::vector<const GenParticle*>
getDaughters(const GenParticle* mother,
             const std::vector<const GenParticle*>& p_vec)
{
  std::vector<const GenParticle*> output{};
  if (mother->D1 > 0) {
    output.reserve(mother->D2 - mother->D1 + 1);
    for (auto idx = mother->D1; idx <= mother->D2; ++idx) {
      const GenParticle* daughter = p_vec.at(idx);
      output.push_back(daughter);
    }
  }
  return output;
}

const GenParticle*
getHardProcessMotherCopy(const GenParticle* p,
                         const std::vector<const GenParticle*>& p_vec)
{
  // is particle itself is hard process particle
  if (isHardProcess(p)) {
    return p;
  }

  // check if any other copies are hard process particles
  const GenParticle* pcopy = p;
  std::unordered_set<const GenParticle*> dup_check;
  while (getPreviousCopy(pcopy, p_vec)) {
    dup_check.insert(pcopy);
    pcopy = getPreviousCopy(pcopy, p_vec);
    if (isHardProcess(pcopy)) {
      return pcopy;
    }
    if (dup_check.count(pcopy)) {
      break;
    }
  }
  return nullptr;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L511-L523
const GenParticle*
getPreviousCopy(const GenParticle* p,
                const std::vector<const GenParticle*>& p_vec)
{
  if (p->M1 == -1) {
    throw std::runtime_error("got an unexpected M1=-1");
  }

  if (p->M2 != -1) {
    throw std::runtime_error("got an unexpected M2=" + std::to_string(p->M2));
  }

  auto mother = p_vec.at(p->M1);

  return mother->PID == p->PID ? mother : nullptr;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L525-L537
const GenParticle*
getNextCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_vec)
{
  for (const auto& dau : getDaughters(p, p_vec)) {
    if (dau->PID == p->PID) {
      return dau;
    }
  }
  return nullptr;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L392-L404
/////////////////////////////////////////////////////////////////////////////
const GenParticle*
getFirstCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_vec)
{
  const GenParticle* p_copy = p;

  std::unordered_set<const GenParticle*> dup_check;
  // FIXME
  while (getPreviousCopy(p_copy, p_vec) != nullptr) {
    dup_check.insert(p_copy);
    p_copy = getPreviousCopy(p_copy, p_vec);
    if (dup_check.count(p_copy)) {
      return nullptr;
    }
  }
  return p_copy;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L406-L418
const GenParticle*
getLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_vec)
{
  const GenParticle* p_copy = p;
  std::unordered_set<const GenParticle*> dup_check;
  while (getNextCopy(p_copy, p_vec) != nullptr) {
    dup_check.insert(p_copy);
    p_copy = getNextCopy(p_copy, p_vec);
    if (dup_check.count(p_copy)) {
      return nullptr;
    }
  }
  return p_copy;
}

/////////////////////////////////////////////////////////////////////////////
bool
isFirstCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_vec)
{
  return p == getFirstCopy(p, p_vec);
}

bool
isLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_vec)
{
  return p == getLastCopy(p, p_vec);
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L265-L295
// TODO
bool
isHardProcess(const GenParticle* p)
{
  // hard process codes for pythia8 are 21-29 inclusive (currently 21,22,23,24
  // are used)
  return ((p->Status > 20) and (p->Status < 30));
}

bool
isFromHardProcess(const GenParticle* p,
                  const std::vector<const GenParticle*>& p_vec)
{
  return getHardProcessMotherCopy(p, p_vec) != nullptr;
}

const GenParticle*
findDaughter(const GenParticle* mother,
             const int pid,
             const std::vector<const GenParticle*>& p_vec)
{
  mother = getLastCopy(mother, p_vec);
  for (const auto daughter : getDaughters(mother, p_vec)) {
    if (daughter->PID == pid) {
      return daughter;
    }
  }

  return nullptr;
}

const GenParticle*
findLastDaughter(const GenParticle* mother,
                 const int pid,
                 const std::vector<const GenParticle*>& p_vec)
{
  const GenParticle* daughter = findDaughter(mother, pid, p_vec);
  if (daughter) {
    daughter = getLastCopy(daughter, p_vec);
  }
  return daughter;
}

// match parton to jet with shortest distance
// starting with the shortest distance available
// apply some outlier rejection if desired
// https://github.com/cms-sw/cmssw/blob/CMSSW_14_1_0_pre4/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L149-L196
std::vector<int32_t>
matchPartonToJet(const std::vector<const GenParticle*> parton_vec,
                 const std::vector<const Jet*> jet_vec,
                 const float max_distance)
{
  // FIXME update comments
  // prepare vector of pairs with distances between
  // all partons to all jets in the input vectors
  std::vector<std::tuple<double, size_t, size_t>> match_vec{};
  match_vec.reserve(parton_vec.size() * jet_vec.size());

  for (size_t parton_idx = 0; parton_idx < parton_vec.size(); ++parton_idx) {
    for (size_t jet_idx = 0; jet_idx < jet_vec.size(); ++jet_idx) {

      const double distance = VectorUtil::DeltaR(
        jet_vec.at(jet_idx)->P4(), parton_vec.at(parton_idx)->P4());
      match_vec.emplace_back(distance, parton_idx, jet_idx);
    }
  }
  std::sort(match_vec.begin(), match_vec.end());

  std::set<size_t> seen{};

  std::vector<int32_t> parton_jet_idx(parton_vec.size(), -1);

  while (not match_vec.empty()) {
    auto [distance, parton_idx, jet_idx] = match_vec.at(0);
    // use primitive outlier rejection if desired
    if (distance > max_distance) {
      jet_idx = -1;
    }

    if (seen.contains(parton_idx)) {
      throw std::runtime_error("got a parton already seen");
    }

    parton_jet_idx.at(parton_idx) = jet_idx;
    seen.insert(parton_idx);

    // remove all values for the matched parton
    // and the matched jet
    size_t other_parton_idx;
    size_t other_jet_idx;
    // XXX
    // the original code use unsigned int
    for (size_t other_idx = 0; other_idx < match_vec.size(); ++other_idx) {
      std::tie(std::ignore, other_parton_idx, other_jet_idx) =
        match_vec.at(other_idx);
      if ((other_parton_idx == parton_idx) or (other_jet_idx == jet_idx)) {
        match_vec.erase(match_vec.begin() + other_idx);
        --other_idx;
      }
    }
  }

  return parton_jet_idx;
}

// --------------------------------------------------------------------------
// Function: IsDescendantOf
// Checks if 'particle' comes from 'ancestor_candidate' using a const vector.
// --------------------------------------------------------------------------
bool
isDescendantOf(const GenParticle* particle,
               const GenParticle* ancestor_candidate,
               const std::vector<const GenParticle*>& particle_vec)
{

  if (particle == nullptr) {
    throw std::runtime_error("particle is nullptr");
  }

  if (ancestor_candidate == nullptr) {
    throw std::runtime_error("ancestor_candidate is nullptr");
  }

  if (particle == ancestor_candidate) {
    spdlog::warn(
      "IsDescendantOf called with identical particle and ancestor_candidate");
    return false;
  }

  // 2. Setup Search Queue
  std::vector<int> mothers_to_check{};
  std::set<int> seen{};

  // Load immediate mothers
  if (particle->M1 >= 0) {
    int m2 = (particle->M2 >= 0) ? particle->M2 : particle->M1;
    for (int i = particle->M1; i <= m2; ++i) {
      mothers_to_check.push_back(i);
    }
  }

  // 3. Loop Upwards
  while (not mothers_to_check.empty()) {
    const int idx = mothers_to_check.back();
    mothers_to_check.pop_back();

    // Safety: Check bounds and cycles
    if ((idx < 0) or (idx >= (int)particle_vec.size())) {
      continue;
    }

    if (seen.count(idx)) {
      throw std::runtime_error("cycle detected in ancestry chain");
    }
    seen.insert(idx);

    // Access the mother from the vector
    const auto mother = particle_vec.at(idx);

    // --- CHECK ---
    if (mother == ancestor_candidate) {
      return true;
    }

    // Add parents to stack
    if (mother->M1 >= 0) {
      const int parentM2 = (mother->M2 >= 0) ? mother->M2 : mother->M1;
      for (int i = mother->M1; i <= parentM2; ++i) {
        mothers_to_check.push_back(i);
      }
    }
  }

  return false;
}
