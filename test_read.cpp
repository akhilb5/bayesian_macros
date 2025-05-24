#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
const double epsilon = 1e-10;

std::vector<std::pair<double, double>> Data;
std::vector<std::pair<double, double>> GSResponse;
std::vector<std::pair<double, double>> Level1a;
std::vector<std::pair<double, double>> Level1b;

std::vector<std::pair<double, double>> read_tsv(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<double, double>> data;

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    std::string line;
//    int line_count = 0;
    while (std::getline(file, line)) {
//        if (line_count >= 8192) break;
        std::istringstream iss(line);
        double x;
        double y;
        if (iss >> x >> y) {
            data.push_back(std::make_pair(x, y));
        }
//        ++line_count;
    }

    return data;
}

// Simple function to load the data
void test_read() {
    Data = read_tsv("Data.txt");
    GSResponse = read_tsv("GSResponse.txt");
    Level1a = read_tsv("Level1a.txt");
    Level1b = read_tsv("Level1b.txt");

    std::cout << "Loaded Data" << Data.size() << " entries.\n";
    std::cout << "Loaded GSResponse" << GSResponse.size() << " entries.\n";
    std::cout << "Loaded Level1a" << Level1a.size() << " entries.\n";
    std::cout << "Loaded Level1b" << Level1b.size() << " entries.\n";
    std::vector<double> d;
    std::vector<double> R0;
    std::vector<double> R1;
    std::vector<double> R2;
    for(int j=0;j<8192; j++) {
        d.push_back(Data[j].second );
        R0.push_back(GSResponse[j].second );
        R1.push_back(Level1a[j].second );
        R2.push_back(Level1b[j].second );
    }
    // 1. Create a histogram
    auto *hData = new TH1D("hData", "Data",8192, 0, 8192);
    auto *hGSResponse = new TH1D("hGSResponse", "GSResponse",8192, 0, 8192);
    auto *hLevel1a = new TH1D("hLevel1a", "Level1a",8192, 0, 8192);
    auto *hLevel1b = new TH1D("hLevel1b", "Level1b",8192, 0, 8192);
    //auto *added = new TH1D("added", "added",8192, 0, 8192);
    auto *scaled_hGsResponse = new TH1D("scled_hGSResponse", "scled_hGSResponse",8192, 0, 8192);
    auto *scaled_hLevel1a = new TH1D("scled_hLevel1a", "scled_hLevel1a",8192, 0, 8192);
    auto *scaled_hLevel1b = new TH1D("scled_hLevel1b", "scled_hLevel1b",8192, 0, 8192);

    // 2. Fill the histogram
    for (size_t i = 0; i < d.size(); ++i) {
        hData->Fill(i,d[i] );
    }
    for(size_t i = 0; i < R0.size(); ++i) {
        hGSResponse->Fill(i,R0[i] );
        hLevel1a->Fill(i,R1[i]);
        hLevel1b->Fill(i,R2[i]);
    }
    auto *c1 = new TCanvas("c1", "Data and responses ", 800, 600);
    // 3. Draw it
    hData->Draw("HIST");
    hGSResponse->Draw("same:HIST");
    hLevel1a->Draw("same:HIST");
    hLevel1b->Draw("same:HIST");
    
    std::vector<double> s0_history;
    std::vector<double> s1_history;
    std::vector<double> s2_history;


    double s[3] = {0.33,0.33,0.33};
    double sumRa[3] = {0.0, 0.0, 0.0};
    for (int k=360; k<1054; k++) {
        sumRa[0] += R0[k];
        sumRa[1] += R1[k];
        sumRa[2] += R2[k];
    }
    //std::cout << "sumRa[0] = " << sumRa[0] << ", sumRa[1] = " << sumRa[1] << ", sumRa[2] = " << sumRa[2] << std::endl;
    std::cout << "s[0] = " << s[0]<< ", s[1] = " << s[1]<< ", s[2] = " << s[2] << std::endl;
    for(int r=0; r<100; r++) {
        double sumRsd[3] = {0.0, 0.0, 0.0};
        double sumRs = 0.0;
        for (int i=360; i<1054; i++){
            sumRs = (R0[i]*s[0]) + (R1[i]*s[1]) + (R2[i]*s[2]);
            if (std::abs(sumRs) < epsilon) continue;
            sumRsd[0] += (R0[i] * s[0] * d[i])/sumRs;
            sumRsd[1] += (R1[i] * s[1] * d[i])/sumRs;
            sumRsd[2] += (R2[i] * s[2] * d[i])/sumRs;
            //std::cout << "sumRsd[0] = " << sumRsd[0] << ", sumRsd[1] = " << sumRsd[1] << ", sumRsd[2] = " << sumRsd[2] << std::endl;
        }
        s[0] = sumRsd[0] / sumRa[0];
        s[1] = sumRsd[1] / sumRa[1];
        s[2] = sumRsd[2] / sumRa[2];
        double Sum = s[0] + s[1] + s[2];
        s0_history.push_back(s[0]/Sum);
        s1_history.push_back(s[1]/Sum);
        s2_history.push_back(s[2]/Sum);

        //std::cout << "sumRsd[0] = " << sumRsd[0] << ", sumRsd[1] = " << sumRsd[1] << ", sumRsd[2] = " << sumRsd[2] << ", sumRs = "<<sumRs << std::endl;
        std::cout << "s[0] = " << s[0]/Sum << ", s[1] = " << s[1]/Sum << ", s[2] = " << s[2]/Sum << std::endl;
    }
    std::cout << "s[0] = " << s[0]<< ", s[1] = " << s[1]<< ", s[2] = " << s[2]<< std::endl;
    auto *c2 = new TCanvas("c2", "data,scaled responses and sum", 800, 600);
    // // double st[3] = {0.053, 0.00058, 0.947};
    // double st[3] = {1.0, 1.0, 1.0};
    // for (int i=0; i<3; i++) {
    //     s[i] = st[i];
    // }
    double Sum2 = 1.0;
    //double Sum2 = s[0] + s[1] + s[2];
    TH1D *added = (TH1D*)hGSResponse->Clone("added");
    added->Scale(s[0]/Sum2);
    added->Add(hLevel1a,s[1]/Sum2);
    added->Add(hLevel1b,s[2]/Sum2);
    added->SetLineColor(kRed);
    added->SetLineStyle(2);
    scaled_hGsResponse->Add(hGSResponse,hGSResponse,s[0]/Sum2,-1);
    scaled_hLevel1a->Add(hLevel1a,hLevel1a,s[1]/Sum2,-1);
    scaled_hLevel1b->Add(hLevel1b,hLevel1b,s[2]/Sum2,-1);
    scaled_hGsResponse->SetLineColor(kBlue);
    scaled_hLevel1a->SetLineColor(kGreen+2);
    scaled_hLevel1b->SetLineColor(kOrange+2);
    hData->Draw("HIST");
    scaled_hGsResponse->Draw("same:HIST");
    //hGSResponse->Draw("same:HIST");
    //hLevel1a->Draw("same:HIST");
    //hLevel1b->Draw("same:HIST");
    scaled_hLevel1a->Draw("same:HIST");
    scaled_hLevel1b->Draw("same:HIST");
    added->Draw("same:HIST");
    std::cout << "Integrals" << "Data" << hData->Integral() << ", GSResponse" << hGSResponse->Integral() << ", Level1a" << hLevel1a->Integral() << ", Level1b" << hLevel1b->Integral() << std::endl;
    std::cout << "Integrals" << "scaled_hGSResponse" << scaled_hGsResponse->Integral() << ", scaled_hLevel1a" << scaled_hLevel1a->Integral() << ", scaled_hLevel1b" << scaled_hLevel1b->Integral() << std::endl;    
    std::cout << "Integrals" << "added" << added->Integral() << std::endl;
    auto *g_s0 = new TGraph(s0_history.size());
    auto *g_s1 = new TGraph(s1_history.size());
    auto *g_s2 = new TGraph(s2_history.size());

    for (int i = 0; i < s0_history.size(); ++i) {
        g_s0->SetPoint(i, i, s0_history[i]);
        g_s1->SetPoint(i, i, s1_history[i]);
        g_s2->SetPoint(i, i, s2_history[i]);
    }

    // Style the graphs a little
    g_s0->SetLineColor(kRed);
    g_s1->SetLineColor(kBlue);
    g_s2->SetLineColor(kGreen+2);

    g_s0->SetLineWidth(2);
    g_s1->SetLineWidth(2);
    g_s2->SetLineWidth(2);
    
    g_s0->GetYaxis()->SetRangeUser(0.0, 1.0);

    // Create a canvas
    auto *c3 = new TCanvas("c3", "Scaling Factors vs Iterations", 800, 600);
    // Draw
    g_s0->Draw("AL");  // A = axis, L = line
    g_s1->Draw("L SAME");
    g_s2->Draw("L SAME");

    // Add legend
    auto *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(g_s0, "s0", "l");
    legend->AddEntry(g_s1, "s1", "l");
    legend->AddEntry(g_s2, "s2", "l");
    legend->Draw();

            
}

