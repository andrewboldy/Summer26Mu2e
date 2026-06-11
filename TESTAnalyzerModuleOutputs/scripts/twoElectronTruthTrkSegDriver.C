//----------------------------------------------------------------------------------
//
// twoElectronTruthTrkSegDriver.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   First-stage driver for the modular twoElectronTruthTrkSegAnalyzer workflow.
//   For now, this macro only runs the MCTruthPrinter art module through the
//   twoElectronTruthTrkSegAnalyzer.fcl configuration.
//
//   The intent is to keep this wrapper thin:
//     1. point at the configured FHiCL file,
//     2. supply the input art file or filelist,
//     3. launch mu2e, and
//     4. leave room for later stages that will add more modules.
//
//   In other words, this file is a launcher, not the analysis itself.  The
//   actual event logic still lives in MCTruthPrinter_module.cc.
//
// Usage from ROOT:
//   .L CreatedCode/AnalyzerScripts/twoElectronTruthTrkSegDriver.C+
//   twoElectronTruthTrkSegDriver("path/to/input.root")
//
// Optional arguments:
//   twoElectronTruthTrkSegDriver("path/to/input.root", 100)
//   twoElectronTruthTrkSegDriver("path/to/input.root", -1,
//                                  "Offline/BoldyAnalyzers/fcl/twoElectronTruthTrkSegAnalyzer.fcl")
//
// Notes:
//   - This macro assumes the MCTruthPrinter plugin has been built and is
//     available in the Mu2e art environment.
//   - The default FHiCL file is the one copied into Offline/BoldyAnalyzers/fcl.
//   - The FHiCL currently prints truth particles from the configured
//     SimParticleCollection and uses the selection switches in the file.
//   - nEvents is optional.  A negative value means "let art run the full input".
//   - fclPath can be overridden later when we add more staged configurations.
//
// Example command that runs the driver:
//   root -l -q "CreatedCode/AnalyzerScripts/twoElectronTruthTrkSegDriver.C++(\"/exp/mu2e/data/users/aboldy/Summer26/MDC2025Sim/B2BCeEndpointDiffEnergy/HundredThousandEvents/nts.owner.description.version.sequencer.root\",-1,\"$BOLDY_WORK_DIR/Offline/BoldyAnalyzers/fcl/twoElectronTruthTrkSegAnalyzer.fcl\")"
//
//----------------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <string>

#include <TSystem.h>

namespace
{
  // ROOT input can be either one art file or a text filelist.  The Mu2e
  // launcher uses -s for a single art file and -S for a filelist, so the driver
  // needs to distinguish them.
  bool hasRootSuffix(const std::string& path)
  {
    const std::string suffix = ".root";

    if (path.size() < suffix.size())
    {
      return false;
    }

    return path.substr(path.size() - suffix.size()) == suffix;
  }

  // Build the mu2e command line in one place so the macro stays easy to extend.
  // Keeping the command assembly here means later stages can reuse the same
  // launcher pattern and only add more options when needed.
  std::string buildMu2eCommand(const std::string& inputName,
                               int nEvents,
                               const std::string& fclPath)
  {
    std::ostringstream cmd;

    // art job invocation:
    //   -c selects the FHiCL configuration
    //   -s points at a single art ROOT file
    //   -S points at a filelist of ROOT files
    //   -n is only added when the caller asks for a finite event limit
    cmd << "mu2e -c " << fclPath
        << (hasRootSuffix(inputName) ? " -s " : " -S ")
        << inputName;

    if (nEvents >= 0)
    {
      cmd << " -n " << nEvents;
    }

    return cmd.str();
  }
}

int twoElectronTruthTrkSegDriver(const std::string& inputName,
                                 int nEvents = -1,
                                 const std::string& fclPath = "Offline/BoldyAnalyzers/fcl/twoElectronTruthTrkSegAnalyzer.fcl")
{
  // Basic input validation.
  // A missing input file should fail early with a clear message instead of
  // launching art and failing later with a less useful error.
  if (inputName.empty())
  {
    std::cerr << "twoElectronTruthTrkSegDriver: input file or filelist is empty\n";
    return 1;
  }

  // Confirm the FHiCL file is reachable before we build the job command.
  // The macro defaults to the package-local configuration, but the path is
  // still overridable for later staged workflows.
  if (gSystem->AccessPathName(fclPath.c_str()))
  {
    std::cerr << "twoElectronTruthTrkSegDriver: cannot find FHiCL file: "
              << fclPath << '\n';
    return 2;
  }

  // First-stage execution: run only MCTruthPrinter using the configured FHiCL.
  // Later stages can add more modules, but the launch pattern stays the same.
  // This keeps the driver stable while the analysis chain grows underneath it.
  const std::string command = buildMu2eCommand(inputName, nEvents, fclPath);

  // Emit the final command line so a user can reproduce the exact job launch.
  std::cout << "twoElectronTruthTrkSegDriver: running\n  " << command << '\n';

  // Hand execution over to the Mu2e environment.
  // gSystem->Exec returns the shell exit code from the command it ran.
  const int status = gSystem->Exec(command.c_str());

  if (status != 0)
  {
    // Nonzero status means mu2e failed or the shell command could not be
    // executed successfully.
    std::cerr << "twoElectronTruthTrkSegDriver: mu2e returned status "
              << status << '\n';
  }

  return status;
}

