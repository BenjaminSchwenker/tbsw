#include <fstream>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH2F.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;

int main(int, char **)
// void DrawBoxes()
{

  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // display mode
  gStyle->SetPalette(1);
  gStyle->SetOptStat(11111111);

  // Read config file
  //------------------
  TEnv *mEnv = new TEnv("x0.cfg");

  // u minimum and v maximum values (in mm)
  std::vector<double> umin;
  std::vector<double> vmin;
  std::vector<double> umax;
  std::vector<double> vmax;

  // Determine the number of fit functions from the cfg file
  int num_MA = 0; // Number of additional measurement areas
  int num_line =
      0; // Number of lines, which are used for BE gradient calibration
  int num_rectangle = 0; // Number of rectangular areasle, which include a whole
                         // grid of measurement areas

  TString MAname;
  MAname.Form("MA%i", num_MA + 1);

  TString linename;
  linename.Form("line%i", num_line + 1);

  TString rectanglename;
  rectanglename.Form("rectangle%i", num_rectangle + 1);

  while (mEnv->GetValue(MAname + ".exist", 0) != 0) {
    umin.push_back(mEnv->GetValue(MAname + ".ucenter", 0.0) -
                   0.5 * mEnv->GetValue(MAname + ".ulength", 0.0));
    umax.push_back(mEnv->GetValue(MAname + ".ucenter", 0.0) +
                   0.5 * mEnv->GetValue(MAname + ".ulength", 0.0));
    vmin.push_back(mEnv->GetValue(MAname + ".vcenter", 0.0) -
                   0.5 * mEnv->GetValue(MAname + ".vlength", 0.0));
    vmax.push_back(mEnv->GetValue(MAname + ".vcenter", 0.0) +
                   0.5 * mEnv->GetValue(MAname + ".vlength", 0.0));

    num_MA++;
    MAname.Form("MA%i", num_MA + 1);
  }

  while (mEnv->GetValue(linename + ".exist", 0) != 0) {
    // Find number of measurement areas along line, typically one steplength
    // parameter should be 0 and the other non-zero
    if (mEnv->GetValue(linename + ".usteplength", 0.0) != 0.0) {

      for (double d = mEnv->GetValue(linename + ".startu", 0.0);
           d < (mEnv->GetValue(linename + ".startu", 0.0) +
                mEnv->GetValue(linename + ".ulength", -1.0));
           d += mEnv->GetValue(linename + ".usteplength", 10.0)) {
        // Compute/Read out measurement area parameters
        umin.push_back(d -
                       0.5 * mEnv->GetValue(linename + ".usteplength", 0.0));
        umax.push_back(d +
                       0.5 * mEnv->GetValue(linename + ".usteplength", 0.0));
        vmin.push_back(mEnv->GetValue(linename + ".startv", 0.0) -
                       0.5 * mEnv->GetValue(linename + ".vlength", 0.0));
        vmax.push_back(mEnv->GetValue(linename + ".startv", 0.0) +
                       0.5 * mEnv->GetValue(linename + ".vlength", 0.0));
      }

    }

    else if (mEnv->GetValue(linename + ".vsteplength", 0.0) != 0.0) {
      for (double d = mEnv->GetValue(linename + ".startv", 0.0);
           d < (mEnv->GetValue(linename + ".startv", 0.0) +
                mEnv->GetValue(linename + ".vlength", -1.0));
           d += mEnv->GetValue(linename + ".vsteplength", 10.0)) {
        // Compute/Read out measurement area parameters
        umin.push_back(mEnv->GetValue(linename + ".startu", 0.0) -
                       0.5 * mEnv->GetValue(linename + ".ulength", 0.0));
        umax.push_back(mEnv->GetValue(linename + ".startu", 0.0) +
                       0.5 * mEnv->GetValue(linename + ".ulength", 0.0));
        vmin.push_back(d -
                       0.5 * mEnv->GetValue(linename + ".vsteplength", 0.0));
        vmax.push_back(d +
                       0.5 * mEnv->GetValue(linename + ".vsteplength", 0.0));
      }
    }

    else
      continue;

    num_line++;
    linename.Form("line%i", num_line + 1);
  }

  while (mEnv->GetValue(rectanglename + ".exist", 0) != 0) {
    // Find number of measurement areas in rectangle
    if ((mEnv->GetValue(rectanglename + ".usteplength", 0.0) != 0.0) &&
        (mEnv->GetValue(rectanglename + ".usteplength", 0.0) != 0.0)) {

      for (double d_u = mEnv->GetValue(rectanglename + ".startu", 0.0);
           d_u < (mEnv->GetValue(rectanglename + ".startu", 0.0) +
                  mEnv->GetValue(rectanglename + ".ulength", -1.0));
           d_u += mEnv->GetValue(rectanglename + ".usteplength", 10.0)) {
        for (double d_v = mEnv->GetValue(rectanglename + ".startv", 0.0);
             d_v < (mEnv->GetValue(rectanglename + ".startv", 0.0) +
                    mEnv->GetValue(rectanglename + ".vlength", -1.0));
             d_v += mEnv->GetValue(rectanglename + ".vsteplength", 10.0)) {
          // Compute/Read out measurement area parameters
          umin.push_back(
              d_u - 0.5 * mEnv->GetValue(rectanglename + ".usteplength", 0.0));
          umax.push_back(
              d_u + 0.5 * mEnv->GetValue(rectanglename + ".usteplength", 0.0));
          vmin.push_back(
              d_v - 0.5 * mEnv->GetValue(rectanglename + ".vsteplength", 0.0));
          vmax.push_back(
              d_v + 0.5 * mEnv->GetValue(rectanglename + ".vsteplength", 0.0));
        }
      }

    } else
      continue;

    num_rectangle++;
    rectanglename.Form("rectangle%i", num_rectangle + 1);
  }

  // TString for the input root file name
  TString filename, histoname, range;
  filename = "X0image";

  // Copy the X0 Analysis Root file
  TFile *rootFile = new TFile(filename, "READ");

  rootFile->cd("");
  gStyle->SetPalette(1, 0);

  // X0 map
  TH2F *hX0map = (TH2F *)rootFile->Get("x0_image");

  // Create a new canvas
  TCanvas *c = new TCanvas("c", "c", 1400, 1000);
  hX0map->SetMaximum(2.3);
  hX0map->SetStats(kFALSE);
  hX0map->Draw("colz");

  hX0map->GetZaxis()->SetTitle("X/X_{0}[%]");
  hX0map->GetZaxis()->SetLabelSize(0.07f);
  hX0map->GetZaxis()->SetTitleSize(0.07f);
  hX0map->GetZaxis()->SetTitleOffset(0.10f);

  hX0map->GetXaxis()->SetTitle("u[mm]");
  hX0map->GetXaxis()->SetLabelSize(0.07f);
  hX0map->GetXaxis()->SetTitleSize(0.07f);
  hX0map->GetXaxis()->SetTitleOffset(0.6f);

  hX0map->GetYaxis()->SetTitle("v[mm]");
  hX0map->GetYaxis()->SetLabelSize(0.07f);
  hX0map->GetYaxis()->SetTitleSize(0.07f);
  hX0map->GetYaxis()->SetTitleOffset(0.5);

  hX0map->SetTitle("");

  for (size_t i = 0; i < umin.size(); i++) {
    TBox *box = new TBox(umin[i], vmin[i], umax[i], vmax[i]);
    TString number;
    number.Form("%lu", i + 1);
    TText *text = new TText(
        0.5 * (umax[i] + umin[i]) - 0.25 * (umax[i] - umin[i]),
        0.5 * (vmax[i] + vmin[i]) - 0.25 * (vmax[i] - vmin[i]), number);
    text->SetTextColor(1);

    if (i < 9) {
      text->DrawText(0.5 * (umax[i] + umin[i]) - 0.25 * (umax[i] - umin[i]),
                     0.5 * (vmax[i] + vmin[i]) - 0.25 * (vmax[i] - vmin[i]),
                     number);
    }

    else {
      text->DrawText(0.5 * (umax[i] + umin[i]) - 0.42 * (umax[i] - umin[i]),
                     0.5 * (vmax[i] + vmin[i]) - 0.25 * (vmax[i] - vmin[i]),
                     number);
    }

    box->SetFillStyle(0);
    box->SetLineColor(1);
    box->SetLineWidth(2);
    box->Draw();
  }

  c->Modified();
  c->Update();

  TPaletteAxis *palette =
      (TPaletteAxis *)hX0map->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.9025);
  palette->SetX2NDC(0.925);

  c->SaveAs(filename + "_Boxes.pdf");

  return 0;
}
