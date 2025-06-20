#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <cmath>

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
namespace fs = std::filesystem;

const double epsilon = 1e-10;

std::vector<std::pair<double, double>> read_tsv(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<double, double>> data;
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y) {
            data.push_back(std::make_pair(x, y));
        }
    }
    return data;
}

void process_directory(int fit_min, int fit_max, const std::string& directory) {
    std::vector<std::pair<double, double>> Data;
    std::vector<std::vector<double>> ResponseData;
    std::vector<std::string> responseFilenames;

    // 1. Read all files in directory
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".txt") {
            std::string filename = entry.path().filename().string();
            if (filename == "Data.txt") {
                Data = read_tsv(entry.path().string());
                std::cout << "Loaded Data: " << Data.size() << " entries.\n";
            } else {
                std::vector<std::pair<double, double>> temp = read_tsv(entry.path().string());
                std::vector<double> values;
                for (const auto& p : temp)
                    values.push_back(p.second);
                ResponseData.push_back(values);
                responseFilenames.push_back(filename);
                std::cout << "File: " << filename << " has " << values.size() << " entries.\n";
            }
        }
    }
    std::ofstream txt("../response_file_names_order.txt");
    for (int i = 0; i < responseFilenames.size(); ++i) {
        std::string filename = responseFilenames[i];
        txt << i << "\t" << filename << "\n";
    }
    size_t N = ResponseData[0].size();//Data.size()/2;
    size_t R = ResponseData.size();
    std::vector<double> d;
    for (size_t i = 0; i < N; ++i)
        d.push_back(Data[i].second);

    // 2. Build Data Histogram
    auto* hData = new TH1D("hData", "Data", N, 0, N);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(2);
    for (size_t i = 0; i < N; ++i)
        hData->SetBinContent(i + 1, d[i]);

    // 3. Compute sumRa
    std::vector<double> sumRa(R, 0.0);
    //for (size_t k = 360; k < 1054; ++k)
    for (size_t k = fit_min; k < fit_max; ++k)
        for (size_t i = 0; i < R; ++i)
            sumRa[i] += ResponseData[i][k];

    // 4. EM Algorithm
    std::vector<double> s(R, 1);
    std::vector<std::vector<double>> s_history(R);
    for (int iter = 0; iter < 100; ++iter) {
        std::vector<double> sumRsd(R, 0.0); //15,640
        //for (int i = 360; i < 1054; ++i) {
        for (int i = fit_min; i < fit_max; ++i) {
            //650 T 620 C
            double sumRs = 0.0;
            for (size_t j = 0; j < R; ++j)
                sumRs += ResponseData[j][i] * s[j];
            //std::cout << "sumRs = " << sumRs << "sumRsd[7]"<< sumRsd[7]<< "\n";
            if (std::abs(sumRs) < epsilon) continue;
            for (size_t j = 0; j < R; ++j)
                sumRsd[j] += ResponseData[j][i] * s[j] * d[i] / sumRs;
        }

        //double sum_s = 0.0;
        for (size_t j = 0; j < R; ++j) {
            s[j] = sumRsd[j] / sumRa[j];
            //sum_s += s[j];
        }

        for (size_t j = 0; j < R; ++j) {
        //     s[j] /= sum_s;
            s_history[j].push_back(s[j]);
        }
    }
    double sum_s = 0.0;
    for (size_t j = 0; j < R; ++j) {
        sum_s += s[j];
    }

    std::cout << "========== Final Scale Factors ==========\n";
    for (size_t j = 0; j < R; ++j) {
        std::cout << "\033[1;34m" << responseFilenames[j] << "\033[0m"  // Blue for filename
        << " \033[1;32mFinal s[" << j << "] = \033[0m"       // Green for "Final s[...] ="
        << "\033[1;36m" << s[j] << "\033[0m" << "\n";     // Cyan for the scaling factor
    }
    std::cout << "//////////////////////////////////////////////"<<"\n"; 
    for (size_t j = 0; j < R; ++j) {
        if(s[j]> 1e-4) std::cout << "\033[1;34m" << responseFilenames[j] << "\033[0m"  // Blue for filename
                  << " \033[1;32mFinal s[" << j << "] = \033[0m"       // Green for "Final s[...] ="
                  << "\033[1;36m" << s[j] << "\033[0m" << "\n";        // Cyan for the scaling factor
    }
    std::cout << "////////////////////sumS = "<< sum_s <<"//////////////////////////"<<"\n"; 
    for (size_t j = 0; j < R; ++j) {
        if(s[j]> 1e-4) std::cout << "\033[1;34m" << responseFilenames[j] << "\033[0m"  // Blue for filename
                  << " \033[1;32mFinal I[" << j << "] = \033[0m"       // Green for "Final s[...] ="
                  << "\033[1;36m" << 100*(s[j]/sum_s)<< "%" << "\033[0m" << "\n";        // Cyan for the scaling factor
    }

    // 5. Create scaled response histograms and combine them
    std::vector<TH1D*> scaledHists;
    TH1D* sumHist = nullptr;
    //TRandom3* randGen = new TRandom3();

    for (size_t j = 0; j < R; ++j) {
        std::string name = "response_" + std::to_string(j);
        auto* hTemp = new TH1D(name.c_str(), name.c_str(), N, 0, N);
        for (size_t bin = 0; bin < N; ++bin)
            hTemp->SetBinContent(bin + 1, ResponseData[j][bin]);

        // auto* scaled = RescaleHistogram(hTemp, s[j], randGen);
        // scaled->SetLineColor(j % 8 + 1);
        // scaled->SetLineWidth(2);
        // scaledHists.push_back(scaled);
        TH1D *scaled = (TH1D*)hTemp->Clone();
        scaled->Add(hTemp,hTemp, s[j],-1);
        scaled->SetLineColor(j % 8 + 1);
        scaled->SetLineWidth(2);
        scaledHists.push_back(scaled);

        //if (s[j] > 1e-4) {
            if (!sumHist){
                //sumHist = (TH1*)scaled->Clone("sumHist");
                sumHist = (TH1D*)hTemp->Clone("sumHist");
                sumHist->Scale(s[0]);
            }
            else{
                //sumHist = AddHistograms(sumHist, scaled, 1.0, 1.0, "sumHist");
                sumHist->Add(hTemp,s[j]);
                //sumHist = AddHistograms(sumHist, hTemp, s[j], 1.0, "sumHist");
            }

        //}
    }

    // 6. Draw: Data vs Sum
    auto* cCompare = new TCanvas("cCompare", "Data vs Sum of Responses", 800, 600);
    cCompare->SetLogy();
    hData->Draw("HIST");
    sumHist->SetLineColor(kRed);
    sumHist->SetLineStyle(1);
    sumHist->Draw("HIST SAME");
    scaledHists[43]->Draw("HIST SAME");
    scaledHists[29]->Draw("HIST SAME");
    scaledHists[36]->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
    leg->AddEntry(hData, "Original Data", "l");
    leg->AddEntry(sumHist, "Sum of Scaled Responses", "l");
    leg->Draw();

    // 7. Draw: All Scaled Histograms
    auto* cAll = new TCanvas("cAll", "All Scaled Response Histograms", 800, 600);
    bool first = true;
    cAll->SetLogy();
    hData->Draw("HIST");
    for (size_t j = 0; j < R; ++j) {
        if (s[j] > 1e-8) {
            if (first) {
                scaledHists[j]->Draw("HIST SAME");
                //std::cout<<"scaledintergral"<<scaledHists[j]->Integral(1,-1)<<"\n";
                first = false;
            } else {
                scaledHists[j]->Draw("HIST SAME");
                //std::cout<<"scaledintergral"<<scaledHists[j]->Integral(1,-1)<<"\n";
            }
        }
    }
    sumHist->Draw("HIST SAME");
    // Plot scaling factor history
    auto* canvas = new TCanvas("canvas", "Scaling Factors", 800, 600);
    std::vector<TGraph*> graphs;
    int colors[] = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange, kViolet, kTeal};

    for (size_t j = 0; j < R; ++j) {
        auto* g = new TGraph(s_history[j].size());
        for (int i = 0; i < s_history[j].size(); ++i)
            g->SetPoint(i, i, s_history[j][i]);
        g->SetLineColor(colors[j % 8]);
        g->SetLineWidth(2);
        graphs.push_back(g);
    }
    double yMax = 0.0;
    for (auto* g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (y > yMax) yMax = y;
        }
    }
    graphs[0]->GetYaxis()->SetRangeUser(0, yMax * 1.1);  // Add 10% margin


    graphs[0]->Draw("AL");
    for (size_t j = 1; j < R; ++j) {
        graphs[j]->Draw("L SAME");
    }    
}
int bayes_multiple_response_sum_compare_no_reScale() {
    process_directory(5, 1024, "/Users/akhil/work_dir/baysean_example_UTK/I136gs_txt_Total");
    //process_directory(2, 1024, "/Users/akhil/work_dir/baysean_example_UTK/I136gs_txt_center");
    //process_directory("/Users/akhil/work_dir/baysean_example_UTK/Cs137");
    //process_directory(2,800,"/Users/akhil/work_dir/baysean_example_UTK/I136m_txt_Total");
    //process_directory(2,800,"/Users/akhil/work_dir/baysean_example_UTK/I136m_txt_Total/I136m_full");
    //process_directory(5,1024,"/Users/akhil/work_dir/baysean_example_UTK/I136m_txt_Center");
    return 0;
}