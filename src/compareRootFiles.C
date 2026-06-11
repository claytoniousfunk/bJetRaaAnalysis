// compareRootFiles.C
// Usage: root -l -b -q 'compareRootFiles.C("file1.root","file2.root")'
// Compares all histograms in two ROOT files bin-by-bin and reports differences.

#include <iostream>
#include <vector>
#include <string>

void compareRootFiles(const char* path1, const char* path2, double tol = 1e-6) {

  TFile *f1 = TFile::Open(path1);
  TFile *f2 = TFile::Open(path2);

  if (!f1 || f1->IsZombie()) { std::cout << "ERROR: cannot open " << path1 << "\n"; return; }
  if (!f2 || f2->IsZombie()) { std::cout << "ERROR: cannot open " << path2 << "\n"; return; }

  TList *keys1 = f1->GetListOfKeys();
  TList *keys2 = f2->GetListOfKeys();

  int n_match = 0, n_diff = 0, n_missing = 0;

  // check every key in file1 against file2
  TIter next(keys1);
  TKey *key;
  while ((key = (TKey*) next())) {

    const char *name = key->GetName();
    const char *cls  = key->GetClassName();

    bool is1D = TString(cls).Contains("TH1");
    bool is2D = TString(cls).Contains("TH2");
    if (!is1D && !is2D) continue;

    TObject *obj2 = f2->Get(name);
    if (!obj2) {
      std::cout << "MISSING in file2: " << name << "\n";
      n_missing++;
      continue;
    }

    if (is1D) {
      TH1 *h1 = (TH1*) f1->Get(name);
      TH1 *h2 = (TH1*) f2->Get(name);

      if (h1->GetNbinsX() != h2->GetNbinsX()) {
        std::cout << "BINNING MISMATCH: " << name
                  << "  (" << h1->GetNbinsX() << " vs " << h2->GetNbinsX() << " bins)\n";
        n_diff++;
        continue;
      }

      bool differs = false;
      double max_diff = 0.;
      int max_bin = -1;
      for (int i = 0; i <= h1->GetNbinsX() + 1; i++) {
        double v1 = h1->GetBinContent(i);
        double v2 = h2->GetBinContent(i);
        double diff = std::abs(v1 - v2);
        if (diff > tol) {
          differs = true;
          if (diff > max_diff) { max_diff = diff; max_bin = i; }
        }
      }

      if (differs) {
        std::cout << "DIFF: " << name
                  << "  max_diff=" << max_diff
                  << " at bin " << max_bin
                  << " (pT~" << h1->GetBinCenter(max_bin) << ")\n";
        n_diff++;
      } else {
        n_match++;
      }

    } else { // TH2
      TH2 *h1 = (TH2*) f1->Get(name);
      TH2 *h2 = (TH2*) f2->Get(name);

      if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
        std::cout << "BINNING MISMATCH: " << name << "\n";
        n_diff++;
        continue;
      }

      bool differs = false;
      double max_diff = 0.;
      int max_bx = -1, max_by = -1;
      for (int ix = 0; ix <= h1->GetNbinsX() + 1; ix++) {
        for (int iy = 0; iy <= h1->GetNbinsY() + 1; iy++) {
          double v1 = h1->GetBinContent(ix, iy);
          double v2 = h2->GetBinContent(ix, iy);
          double diff = std::abs(v1 - v2);
          if (diff > tol) {
            differs = true;
            if (diff > max_diff) { max_diff = diff; max_bx = ix; max_by = iy; }
          }
        }
      }

      if (differs) {
        std::cout << "DIFF: " << name
                  << "  max_diff=" << max_diff
                  << " at bin (" << max_bx << "," << max_by << ")\n";
        n_diff++;
      } else {
        n_match++;
      }
    }
  }

  // check for keys in file2 not present in file1
  TIter next2(keys2);
  while ((key = (TKey*) next2())) {
    const char *name = key->GetName();
    const char *cls  = key->GetClassName();
    bool is1D = TString(cls).Contains("TH1");
    bool is2D = TString(cls).Contains("TH2");
    if (!is1D && !is2D) continue;
    if (!f1->Get(name)) {
      std::cout << "MISSING in file1: " << name << "\n";
      n_missing++;
    }
  }

  std::cout << "\n=== Summary ===\n";
  std::cout << "  Matching:  " << n_match   << "\n";
  std::cout << "  Different: " << n_diff    << "\n";
  std::cout << "  Missing:   " << n_missing << "\n";
  if (n_diff == 0 && n_missing == 0)
    std::cout << "  --> Files are identical.\n";
  else
    std::cout << "  --> Files differ.\n";

  f1->Close();
  f2->Close();
}
