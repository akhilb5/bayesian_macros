import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
epsilon = 1e-10


def read_tsv(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    x, y = float(parts[0]), float(parts[1])
                    data.append((x, y))
                except ValueError:
                    continue
    return data


def process_directory(fit_min, fit_max, directory):
    Data = []
    ResponseData = []
    responseFilenames = []

    for entry in os.scandir(directory):
        if entry.name.endswith(".txt"):
            if entry.name == "Data.txt":
                Data = read_tsv(entry.path)
                print(f"Loaded Data: {len(Data)} entries.")
            else:
                temp = read_tsv(entry.path)
                values = [y for x, y in temp]
                ResponseData.append(values)
                responseFilenames.append(entry.name)

    with open("../response_file_names_order.txt", "w") as txt:
        for i, name in enumerate(responseFilenames):
            txt.write(f"{i}\t{name}\n")

    N = len(ResponseData[0])
    R = len(ResponseData)
    d = [Data[i][1] for i in range(N)]

    hData = ROOT.TH1D("hData", "Data", N, 0, N)
    hData.SetLineColor(ROOT.kBlack)
    hData.SetLineWidth(2)
    for i in range(N):
        hData.SetBinContent(i + 1, d[i])

    sumRa = [0.0] * R
    for k in range(fit_min, fit_max):
        for i in range(R):
            sumRa[i] += ResponseData[i][k]

    s = [0.10] * R
    s_history = [[] for _ in range(R)]
    for _ in range(250):
        sumRsd = [0.0] * R
        for i in range(fit_min, fit_max):
            sumRs = sum(ResponseData[j][i] * s[j] for j in range(R))
            if abs(sumRs) < epsilon:
                continue
            for j in range(R):
                sumRsd[j] += ResponseData[j][i] * s[j] * d[i] / sumRs

        for j in range(R):
            s[j] = sumRsd[j] / sumRa[j]
            s_history[j].append(s[j])

    sum_s = sum(s)

    print("========== Final Scale Factors ==========")
    for j in range(R):
        print(f"\033[1;34m{responseFilenames[j]}\033[0m \033[1;32mFinal s[{j}] = \033[0m\033[1;36m{s[j]}\033[0m")
    print("//////////////////////////////////////////////")
    for j in range(R):
        if s[j] > 1e-4:
            print(f"\033[1;34m{responseFilenames[j]}\033[0m \033[1;32mFinal s[{j}] = \033[0m\033[1;36m{s[j]}\033[0m")
    print(f"////////////////////sumS = {sum_s}//////////////////////////")
    for j in range(R):
        if s[j] > 1e-4:
            print(f"\033[1;34m{responseFilenames[j]}\033[0m \033[1;32mFinal I[{j}] = \033[0m\033[1;36m{100 * s[j] / sum_s}%\033[0m")

    scaledHists = []
    sumHist = None

    for j in range(R):
        hTemp = ROOT.TH1D(f"response_{j}", f"response_{j}", N, 0, N)
        for bin in range(N):
            hTemp.SetBinContent(bin + 1, ResponseData[j][bin])

        scaled = hTemp.Clone()
        scaled.Add(hTemp, hTemp, s[j], -1)
        scaled.SetLineColor(j % 8 + 1)
        scaled.SetLineWidth(2)
        scaledHists.append(scaled)

        if not sumHist:
            sumHist = hTemp.Clone("sumHist")
            sumHist.Scale(s[j])
        else:
            sumHist.Add(hTemp, s[j])

    cCompare = ROOT.TCanvas("cCompare", "Data vs Sum of Responses", 800, 600)
    cCompare.SetLogy()
    hData.Draw("HIST")
    sumHist.SetLineColor(ROOT.kRed)
    sumHist.SetLineStyle(1)
    sumHist.Draw("HIST SAME")
    scaledHists[56].Draw("HIST SAME")

    leg = ROOT.TLegend(0.7, 0.6, 0.9, 0.85)
    leg.AddEntry(hData, "Original Data", "l")
    leg.AddEntry(sumHist, "Sum of Scaled Responses", "l")
    leg.Draw()

    cAll = ROOT.TCanvas("cAll", "All Scaled Response Histograms", 800, 600)
    cAll.SetLogy()
    hData.Draw("HIST")
    for j in range(R):
        if s[j] > 1e-8:
            scaledHists[j].Draw("HIST SAME")
    sumHist.Draw("HIST SAME")

    canvas = ROOT.TCanvas("canvas", "Scaling Factors", 800, 600)
    graphs = []
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kCyan + 2, ROOT.kOrange, ROOT.kViolet, ROOT.kTeal]

    for j in range(R):
        g = ROOT.TGraph(len(s_history[j]))
        for i, val in enumerate(s_history[j]):
            g.SetPoint(i, i, val)
        g.SetLineColor(colors[j % len(colors)])
        g.SetLineWidth(2)
        graphs.append(g)

    yMax = max(max(ys) for ys in s_history)
    graphs[0].GetYaxis().SetRangeUser(0, yMax * 1.1)
    graphs[0].Draw("AL")
    for j in range(1, R):
        graphs[j].Draw("L SAME")


def bayes_multiple_response_sum_compare_no_reScale():
    process_directory(5, 1024, "/Users/akhil/work_dir/baysean_example_UTK/I136gs_txt_Total")
    return 0


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(False)
    bayes_multiple_response_sum_compare_no_reScale()
