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
            data.push_back(std::make_pair(x, y));  // replaced emplace_back with push_back
        }
    }
    return data;
}

void process_directory(const std::string& directory) {
    std::vector<std::pair<double, double>> Data;
    std::vector<std::vector<double>> ResponseData;
    std::vector<std::string> responseFilenames;

    // 1. Read all files in directory
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".txt") {
            std::string filename = entry.path().filename().string();
            if (filename == "Data.txt") {
                Data = read_tsv(entry.path().string());
                std::cout << "Loaded Data" << Data.size() << " entries." << "data[13][1]"<< Data[13].first <<"\n";
            } else {
                std::vector<std::pair<double, double>> temp = read_tsv(entry.path().string());
                std::vector<double> values;
                for (const auto& p : temp) {
                    values.push_back(p.second);
                }
                ResponseData.push_back(values);
                responseFilenames.push_back(filename);
                std::cout << "Loaded Response" << filename << "  " <<ResponseData.size() << " entries."<< "ResponseData[13][1]"<< ResponseData[0][13] <<"\n";
            }
        }
    }

    size_t N = Data.size();
    size_t R = ResponseData.size();
    std::cout << "N = " << N << ", R = " << R << " " << ResponseData[1].size() << "\n";
    std::vector<double> d;
    for (size_t i = 0; i < N; ++i)
        d.push_back(Data[i].second);

    // 2. Build Histograms (optional)
    auto* hData = new TH1D("hData", "Data", N, 0, N);
    hData->SetLineColor(kBlack);
    for (size_t i = 0; i < N; ++i)
        hData->Fill(i, d[i]);
    
    auto* c = new TCanvas("hist", "Hist", 800, 600);    
    hData->Draw("HIST");

    // 3. Compute sumRa
    std::vector<double> sumRa(R, 0.0); //15,640
    for (size_t k = 360; k < 1054; k++) {
        for (size_t i = 0; i < R; ++i) {
            sumRa[i] += ResponseData[i][k];
        }
        //std::cout << "sumRa[" << k << "] [4] = " << sumRa[8] << "\n";
    }
    //for (size_t i = 0; i < R; ++i) {std::cout << sumRa[i] << "\n";}

    // 4. EM iterations
    std::vector<double> s(R, 0.33);  // Initialize equally
    std::vector<std::vector<double>> s_history(R);

    for (int iter = 0; iter < 8000; ++iter) {
        std::vector<double> sumRsd(R, 0.0);
        for (int i = 360; i < 1054; i++) {
            double sumRs = 0.0;
            for (size_t j = 0; j < R; ++j)
                sumRs += ResponseData[j][i] * s[j];
                //std::cout << "sumRs" << sumRs << "\n";
            if (std::abs(sumRs) < epsilon) continue;

            for (size_t j = 0; j < R; ++j) {
                sumRsd[j] += ResponseData[j][i] * s[j] * d[i] / sumRs;
                //if(sumRs ==0) std::cout << "sumRsd" << sumRsd[j] << "\n";
            }
        }
        //std::cout << "sumRsd" << sumRsd[0] << "\n";

        double sum_s = 0.0;
        for (size_t j = 0; j < R; ++j) {
            s[j] = sumRsd[j] / sumRa[j];
            sum_s += s[j];
        }

        // Normalize and record history
        for (size_t j = 0; j < R; ++j) {
            //s[j] /= sum_s;
            s_history[j].push_back(s[j]);
        }

        // // Output per iteration
        // std::cout << "Iteration " << iter << ": ";
        // for (size_t j = 0; j < R; ++j) {
        //     std::cout << "s[" << j << "] = " << s[j] << " ";
        // }
        // std::cout << "\n";
    }
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
    // 5. Plot scaling factor history
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

    graphs[0]->Draw("AL");
    for (size_t j = 1; j < R; ++j) {
        if(s[j]> 1e-4)graphs[j]->Draw("L SAME");
    }

    // auto* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // for (size_t j = 0; j < R; ++j) {
    //     legend->AddEntry(graphs[j], responseFilenames[j].c_str(), "l");
    // }
    // legend->Draw();
}
int bayes_multiple_response() {
    //process_directory("/Users/akhil/work_dir/baysean_example_UTK/I136gs_txt");
    process_directory("/Users/akhil/work_dir/baysean_example_UTK/Cs137");
    return 0;
}
