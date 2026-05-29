//----------------------------------------------------------------------------------
//
// CaloHitter.C
//
// Small ROOT macro wrapper around the header-only CaloHitter helper.
//
// Usage from a ROOT-enabled Mu2e shell:
//   root -l -b -q 'CreatedCode/HistogramMakers/CaloHitter.C("CaloHitter_CalorimeterDisks.pdf")'
//
// The implementation lives in CaloHitter.hh so other macros can include the same
// Calorimeter/CaloCrystal structures without needing to compile a separate library.
//
//----------------------------------------------------------------------------------

// Pull in the CaloCrystal, Calorimeter, geometry-building, drawing, and PDF-saving
// functions from the calohitter namespace.
#include "CaloHitter.hh"

// std::string is used for the output PDF path argument.
#include <string>

// ROOT calls this function when the macro is executed.  It returns the canvas so
// an interactive ROOT session can keep inspecting or modifying the drawing after
// the PDF has been written.
TCanvas* CaloHitter(const std::string& outputPdf = "CaloHitter_CalorimeterDisks.pdf")
{
  // Build the blank two-disk calorimeter geometry, draw it, save it as outputPdf,
  // and print the disk/total crystal counts to the terminal.
  return calohitter::saveCalorimeterPdf(outputPdf);
}

