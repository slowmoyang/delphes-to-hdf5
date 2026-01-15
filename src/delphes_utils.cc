#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

#include "Math/VectorUtil.h"
#include "TSystem.h"

// std
#include <cstdlib>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

void
setupDelphes()
{
  const std::string env_name = "CONDA_PREFIX";
  if (not std::getenv(env_name.c_str())) {
    throw std::runtime_error(env_name + " not defined");
  }

  const fs::path prefix{ std::getenv(env_name.c_str()) };
  if (not fs::exists(prefix)) {
    throw std::runtime_error("Prefix does not exist: " + prefix.string());
  }

  fs::path so_file = prefix / "lib" / "libDelphes.so";
  fs::path include_dir = prefix / "include";

  if (not fs::exists(so_file)) {
    throw std::runtime_error("shared object not found: " + so_file.string());
  }

  if (not fs::exists(include_dir)) {
    throw std::runtime_error("include directory not found: " +
                             include_dir.string());
  }

  gInterpreter->AddIncludePath(include_dir.c_str());
  gSystem->Load(so_file.c_str());
  gInterpreter->Declare("#include \"classes/DelphesClasses.h\"");
}
