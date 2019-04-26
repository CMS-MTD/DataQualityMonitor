import numpy as np
import os, re, shutil
import argparse
import itertools
from prettytable import PrettyTable

import ROOT as rt
from root_numpy import tree2array, tree2rec
from lib.histo_utilities import create_TH1D, create_TH2D, quantile, rootTH2_to_np
from lib.cebefo_style import cebefo_style, Set_2D_colz_graphics

donotdelete = []



def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="input root file, if -N is given XX is replaced with runNumber", nargs='+')
    parser.add_argument("-C", "--config", type=str, default='config/FNAL_TB_1811/VME_SiPM.txt', help="Config file")
    parser.add_argument("-S", "--save_loc", type=str, default='./out_plots/', help="Saving location")
    parser.add_argument("-N", "--runs_interval", default=None, help="Runs to run", nargs='+')

    parser.add_argument("--figformat", type=str, default='.png', help="Format of output figures")

    parser.add_argument("-B", "--batch", default=True, action='store_false', help="Root batch mode")
    parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Activate verbose mode")

    args = parser.parse_args()
    return args


class Config:
    def __init__(self, infile):
        self.raw_conf = open(infile, 'r').readlines()
        self.channel = {}
        self.ch_ordered = []
        self.plots = []
        self.labels = []
        self.xy_center = [0,0, 15, 15]

        for i, l in enumerate(self.raw_conf):
            if l[0] == '#':
                continue
            l = l[0:-1].split()
            if '-->Print' in l[0]:
                self.plots = l[1:]
            elif '-->XYcenter' in l[0]:
                self.xy_center = [float(l[1]), float(l[2]), float(l[3]), float(l[4])]
            elif '-->SpaceBinSize' in l[0]:
                self.space_bin_size = float(l[1])
            elif '-->TracksCleaning' in l[0]:
                self.TracksCleaning = {}
                cuts = ''
                for aux in l[1:]:
                    if aux[:7] == 'CutOnCh':
                        self.TracksCleaning['ch'] = int(aux[7:])
                        print 'Cuts from ch', int(aux[7:]), 'will be used for cleaning tracks'
                    else:
                        cuts += aux + ' && '
                if len(cuts) > 3:
                    self.TracksCleaning['cuts'] = cuts[:-4]
                    print 'Cleaning tracks with:', cuts[:-4]
                else:
                    self.TracksCleaning['cuts'] = ''
            elif '-->TracksConsistency' in l[0]:
                self.TracksConsistency = {}
                self.TracksConsistency['Bad'] = 0
                self.TracksConsistency['Checked'] = 0
                for aux in l[1:]:
                    if aux[:5] == 'RefCh':
                        self.TracksConsistency['RefCh'] = int(aux[5:])
                        # print 'Tracks consistency: ch', int(aux[5:]), 'will be checked'

                    if aux == 'List':
                        self.TracksConsistency['List'] = 1
            elif '-->Var' in l[0]:
                if not hasattr(self, 'Var'):
                    self.Var = {}
                self.Var[l[0][7:-1]] = l[1:]
            elif len(self.labels) == 0:
                self.labels = l[1:]
            else:
                n = int(l[0])
                self.ch_ordered.append(n)
                l = l[1:]
                self.channel[n] = {}
                for k, val in zip(self.labels, l):
                    if 'idx' in k:
                        self.channel[n][k] = int(val)
                    else:
                        self.channel[n][k] = val

class Bauble:
    def __init__(self):
        pass

def define_range_around_peak(h, perc = [0.4, 0.4], Range=[0.0, 99999.0]):
    h_aux = h.Clone('h_aux')
    for i in range(1, h_aux.GetNbinsX()+1):
        if h_aux.GetBinCenter(i) < Range[0] or h_aux.GetBinCenter(i) > Range[1]:
            h_aux.SetBinContent(i, 0.)

    SS = rt.TSpectrum()
    n_pks = SS.Search(h_aux, 2, "", 0.2)

    i_max = 0
    max = -1
    y = SS.GetPositionY()
    x = SS.GetPositionX()
    found = False
    for i in range(n_pks):
        if y[i] > max and x > Range[0] and x < Range[1]:
            max = y[i]
            i_max = i
            found = True

    if found:
        n_pk = h.FindBin(SS.GetPositionX()[i_max])
    else:
        n_pk = h_aux.GetMaximumBin()

    thr = perc[0] * h.GetBinContent(n_pk)
    n_low = n_pk
    while h.GetBinContent(n_low) > thr and n_low > 1:
        n_low -= 1
    x_low = h.GetBinLowEdge(n_low)

    thr = perc[1] * h.GetBinContent(n_pk)
    n_up = n_pk
    while h.GetBinContent(n_up) > thr and n_up < h.GetNbinsX():
        n_up += 1
    x_up = h.GetBinLowEdge(n_up) + h.GetBinWidth(n_up)
    return x_low, x_up, n_pk

def fill_TimeResHisto(k, dt, h2D, out_list, tagin, title_tag, canvas, save_canvas=True):
    q_up, e_up = quantile(1000*dt, 0.15)
    q_dwn, e_dwn = quantile(1000*dt, 0.85)
    if q_up == None or q_dwn == None:
        print '[WARNING] Quantile estimation failed'
        return
    disp_est = 0.5*np.abs(q_up - q_dwn)
    disp_unc = 0.5*np.hypot(e_up, e_dwn)

    out_list.append([disp_est, disp_unc])
    ix = h2D.GetXaxis().FindBin(0.5*(xu+xd))
    iy = h2D.GetYaxis().FindBin(0.5*(yu+yd))
    ig = h2D.GetBin(ix, iy)
    h2D.SetBinContent(ig, disp_est)
    h2D.SetBinError(ig, disp_unc)

    if save_canvas:
        median = np.percentile(dt, 50)
        width = np.abs(np.percentile(dt, 10) - np.percentile(dt, 90))
        tag = tagin+'_'+str(k)+'_'+str(ib)
        h = create_TH1D(dt, 'h_'+tag,
                        title_tag + ' time resolution - x #in [{:.1f}, {:.1f}], y #in [{:.1f}, {:.1f}]'.format(xd, xu, yd, yu),
                        binning = [ None, median-2*width, median+2*width],
                        axis_title = [var_dT + ' [ns]', 'Events'])
        canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
        h.Fit('gaus', 'LQR','', np.percentile(dt, 20), np.percentile(dt, 98))
        h = h.DrawCopy('LE')

        l = rt.TLatex()
        l.SetTextSize(0.04);
        l.DrawLatexNDC(0.15, 0.85, 'Unbiased dispersion: {:.1f} #pm {:.1f}'.format(disp_est, disp_unc))
        canvas['c_'+tag].Update()
        canvas['c_'+tag].SaveAs(out_dir + '/'+tagin+'_bar{:02d}/'.format(k)+tag+'.png')

def Print2D_1D_canvas(l_h2D, l_h1D, outfile, space_lim = None):
    can = rt.TCanvas('c', 'c', 800, 500)
    can.Divide(2,1)
    padL = can.cd(1)
    padL.Divide(1,3)

    pad = padL.cd(1)
    pad.SetRightMargin(0.15)
    l_h2D[0].DrawCopy('colz1')

    pad = padL.cd(2)
    pad.SetRightMargin(0.15)
    l_h2D[1].DrawCopy('colz1')
    pad = padL.cd(3)
    pad.SetRightMargin(0.15)
    l_h2D[2].DrawCopy('colz1')

    pad = can.cd(2)
    leg = rt.TLegend(0.85,0.8,0.95,0.92)

    l_h1D[0].SetLineColor(1)
    l_h1D[0].Draw()
    leg.AddEntry(l_h1D[0], 'L+R', 'l')
    l_h1D[1].SetLineColor(2)
    l_h1D[1].Draw('SAME')
    leg.AddEntry(l_h1D[1], 'L', 'l')
    l_h1D[2].SetLineColor(4)
    l_h1D[2].Draw('SAME')
    leg.AddEntry(l_h1D[2], 'R', 'l')
    leg.Draw()

    can.Update()
    can.SaveAs(outfile)


if __name__ == '__main__':
    cebefo_style()

    rt.gErrorIgnoreLevel = 6000

    args = parsing()

    if args.batch:
        rt.gROOT.SetBatch()

    if args.figformat[0] != '.':
        args.figformat = '.'+args.figformat
    figform = args.figformat

    configurations = Config(args.config)

    if args.runs_interval==None:
        aux = re.search(r'Run[0-9]+', args.input_file[0])
        flag = aux.group(0)[3:]
        if len(args.input_file) > 1:
            aux = re.search(r'Run[0-9]+', args.input_file[-1])
            flag += '_'+aux.group(0)[3:]
    elif len(args.runs_interval)==1 and not args.runs_interval[0].endswith('.txt'):
        N = args.runs_interval[0]
        args.input_file = [args.input_file[0].replace('XXX', str(N))]
        flag = '_'+str(N)
    else:
        if args.runs_interval[0].endswith('.txt'):
            runs = np.sort(np.loadtxt(args.runs_interval[0]).astype(np.int))
            runs = list(runs)
        elif len(args.runs_interval)==2 and args.runs_interval[0] < args.runs_interval[1]:
                runs = range(int(args.runs_interval[0]), int(args.runs_interval[1])+1)
        else:
            runs = args.runs_interval

        file_template = args.input_file[0]
        args.input_file = map(lambda x: file_template.replace('XXX', str(x)), runs)
        if len(runs) > 1:
            flag = '_' + str(runs[0]) + '-' + str(runs[-1])
        else:
            flag = '_' + str(runs[0])

    save_loc = args.save_loc
    if save_loc[-1] != '/':
        save_loc += '/'
    if not os.path.isdir(save_loc):
        os.system('mkdir -p '+save_loc)
        shutil.copy('lib/index.php', save_loc+'index.php')

    out_dir = save_loc + 'Run' + str(flag) + '_' + os.path.basename(args.config[:-4])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        shutil.copy('lib/index.php', out_dir+'/index.php')


    chain = rt.TChain('pulse')
    chain.tree_name = tree_name = 'pulse'
    for i, f in enumerate(args.input_file):
        root_file = rt.TFile.Open( f,"READ");
        if not root_file:
            print "[ERROR]: Input file not found:", f
            continue
        if not (root_file.GetListOfKeys().At(0)):
            print 'Skipping file with no streamer Info!!'
            continue
        tree_name = root_file.GetListOfKeys().At(0).GetName()
        if i == 0:
            if tree_name != 'pulse':
                print 'Change'
                chain = rt.TChain(tree_name)
                chain.tree_name = tree_name
        elif tree_name != chain.tree_name:
            print 'Skipping', f, ' for tree name incompatibility'
            continue
        chain.Add(f)

    canvas = {}
    canvas['amp'] = {}
    canvas['int'] = {}
    canvas['risetime'] = {}
    canvas['wave'] = {}
    canvas['pos'] = {}
    canvas['pos_sel'] = {}
    canvas['amp_w_pos'] = {}
    canvas['int_w_pos'] = {}
    canvas['t_res_raw'] = {}
    canvas['space_corr'] = {}
    canvas['t_res_space'] = {}
    canvas['dt_vs_amp'] = {}
    canvas['dtcorr_vs_amp'] = {}
    canvas['dt_corr'] = {}

    rt.gStyle.SetStatY(0.95)
    rt.gStyle.SetStatX(0.85)
    rt.gStyle.SetStatW(0.15)
    rt.gStyle.SetStatH(0.1)
    rt.gStyle.SetHistLineWidth(2)

    '''=========================== AllTracks ==========================='''
    if 'AllTracks' in configurations.plots:
        name = 'h_alltracks'
        title = 'All tracks at z_DUT[0]'

        var = 'y_dut[0]:x_dut[0]'
        dy = configurations.xy_center[1]
        dx = configurations.xy_center[0]
        width = configurations.xy_center[2]

        h = rt.TH2D(name, title, 100, -width+dx, width+dx, 100, -width+dy, width+dy)
        h.SetXTitle('x [mm]')
        h.SetYTitle('y [mm]')

        sel = ''
        if hasattr(configurations, 'TracksCleaning'):
            sel += configurations.TracksCleaning['cuts']
        chain.Project(name, var, sel)

        canvas['AllTracks'] = rt.TCanvas('c_allTracks', 'c_allTracks', 800, 600)
        h.DrawCopy('colz')

        canvas['AllTracks'].Update()
        canvas['AllTracks'].SaveAs(out_dir + '/AllTracks.png')

    print '\n======================= Channels selection loop ==========================='
    for k in configurations.ch_ordered:
        print '---> Channel', k
        conf = configurations.channel[k]

        line = rt.TLine()
        line.SetLineColor(6)
        line.SetLineWidth(2)
        line.SetLineStyle(10)

        selection = '(amp[{nch}] != 0 && integral[{nch}] != 0)'.format(nch=k)

        '''=========================== Amplitude ==========================='''
        if 'Amp' in configurations.plots:
            name = 'h_amp_'+str(k)
            title = 'Amplitude channel '+str(k)

            amp_aux = np.concatenate(list(tree2array( chain, 'amp['+str(k)+']')))

            h = rt.TH1D(name, title, 80, np.percentile(amp_aux, 1), min(990., np.percentile(amp_aux, 99.99)))
            h.SetXTitle('Peak amplitude [mV]')
            h.SetYTitle('Events / {:.1f} mV'.format(h.GetBinWidth(1)))
            chain.Project(name, 'amp['+str(k)+']')

            Range = [0.0, 9999999.0]
            if 'cut' in conf.keys():
                if conf['cut'][0] in ['a', 'A']:
                    if 'min' in conf.keys():
                        Range[0] = float(conf['min'])
                    if 'max' in conf.keys():
                        Range[1] = float(conf['max'])

            x_low, x_up, n_pk = define_range_around_peak(h, [0.1, 0.1], Range)

            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            h.DrawCopy('E')

            gr = rt.TGraph(1)
            gr.SetPoint(0, h.GetBinCenter(n_pk), h.GetBinContent(n_pk))
            gr.SetMarkerStyle(23)
            gr.SetMarkerColor(2)
            gr.Draw("P")

            if h.GetMaximum() > 5*h.GetBinContent(n_pk):
                canvas['amp'][k].SetLogy()

            if 'cut' in conf.keys():
                if conf['cut'][0] in ['a', 'A']:
                    selection += ' && (amp[{}]>{} && amp[{}]<{})'.format(k, x_low, k, x_up)
                    line.DrawLine(x_low, 0, x_low, h.GetMaximum())
                    line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch{:02d}'.format(k) + figform)


        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            int_aux = -np.concatenate(list(tree2array(chain, 'integral['+str(k)+']', 'integral['+str(k)+'] != 0')))
            h = rt.TH1D(name, title, 100, np.percentile(int_aux, 0.1), np.max(int_aux))
            h.SetXTitle('Integral [pC]')
            h.SetYTitle('Events / {:.1f} pC'.format(h.GetBinWidth(1)))
            chain.Project(name, '-integral['+str(k)+']', '-integral['+str(k)+'] != 0')

            Range = [0.0, 9999999.0]
            if 'cut' in conf.keys():
                if conf['cut'][0] in ['i', 'I']:
                    if 'min' in conf.keys():
                        Range[0] = float(conf['min'])
                    if 'max' in conf.keys():
                        Range[1] = float(conf['max'])

            x_low, x_up, n_pk = define_range_around_peak(h, [0.25, 0.3], Range)

            h.DrawCopy('E')

            gr = rt.TGraph(1)
            gr.SetPoint(0, h.GetBinCenter(n_pk), h.GetBinContent(n_pk))
            gr.SetMarkerStyle(23)
            gr.SetMarkerColor(2)
            gr.Draw("P")

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()

            if 'cut' in conf.keys():
                if conf['cut'][0] in ['i', 'I']:
                    selection += ' && (-integral[{}]>{} && -integral[{}]<{})'.format(k, x_low, k, x_up)
                    line.DrawLine(x_low, 0, x_low, h.GetMaximum())
                    line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch{:02d}'.format(k)+figform)

        '''=========================== Risetime ======================================='''
        if 'Risetime' in configurations.plots:
            print '\tRisetime'
            canvas['risetime'][k] = rt.TCanvas('c_risetime_'+str(k), 'c_int_'+str(k), 800, 600)

            nbins  = 100
            name   = 'h_risetime_' + str(k)
            title  = 'Risetime ch' + str(k)
            rise_aux = (tree2array(chain, 'risetime['+str(k)+']').view(np.recarray).T)[0]
            bin_w = 200.0/1024
            lower_bound = np.percentile(rise_aux, 5) - bin_w*0.5
            upper_bound = np.percentile(rise_aux, 99) + bin_w*0.5
            N_bins = (upper_bound - lower_bound) / bin_w
            h      = rt.TH1D(name, title, int(N_bins), lower_bound, upper_bound)
            h.SetXTitle('Risetime [ns]')
            h.SetYTitle('Events / {:.1f} ns'.format(h.GetBinWidth(1)))
            if conf['idx_ref'] == -1:#do not apply cut again on reference channels
                cut = selection
            else:
                cut = selection + '&&' + configurations.channel[conf['idx_ref']]['sel']

            chain.Project(name, 'risetime['+str(k)+']', cut)

            if (h.GetMaximum() - h.GetMinimum() > 50000):
                canvas['risetime'][k].SetLogy()

            ##ploting histogram
            h.DrawCopy('lE')
            canvas['risetime'][k].Update()
            canvas['risetime'][k].SaveAs(out_dir+'/risetime_ch{:02d}'.format(k)+figform)


        '''=========================== Waveform color chart ==========================='''
        if 'WaveColor' in configurations.plots:
            name = 'h_wave_'+str(k)
            title = 'Waveform color chart channel '+str(k)

            h = rt.TH2D(name, title, 250, 0, 210, 250, -1100, 500)
            h.SetXTitle('Time [ns]')
            h.SetYTitle('Voltage [mV]')

            N_max = chain.GetEntries()
            N_max = max(500, int(N_max*0.1))
            chain.Project(name, 'channel[{}]:time[{}]'.format(k,conf['idx_time'], 'Entry$ < {}'.format(N_max)))

            h.SetStats(0)
            canvas['wave'][k] = rt.TCanvas('c_wave_'+str(k), 'c_wave_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['wave'][k].SetLogz()
            h.DrawCopy('colz')
            canvas['wave'][k].Update()
            canvas['wave'][k].SaveAs(out_dir + '/Waveform_ch{:02d}'.format(k)+figform)

        '''=========================== Track position ==========================='''
        if conf['idx_dut'] >= 0:
            if hasattr(configurations, 'TracksCleaning'):
                selection += ' && ' + configurations.TracksCleaning['cuts']

                if 'ch' in configurations.TracksCleaning.keys():
                    ch = configurations.TracksCleaning['ch']
                    if ch < k:
                        selection += ' && ' + configurations.channel[ch]['sel']

            var = 'y_dut[{}]:x_dut[{}]'.format(conf['idx_dut'], conf['idx_dut'])
            dy = configurations.xy_center[1]
            dx = configurations.xy_center[0]
            width_x = configurations.xy_center[2]
            width_y = configurations.xy_center[3]
            N_bins = [100, 50]
            if 'PosRaw' in configurations.plots:
                name = 'h_pos_'+str(k)
                title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '], for channel {} selected event'.format(k)

                h = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h.SetXTitle('x [mm]')
                h.SetYTitle('y [mm]')
                h.SetZTitle('Events')

                chain.Project(name, var, selection)

                name = 'h_pos1D_'+str(k)
                title = 'Track distribution at z_DUT[' + str(conf['idx_dut']) + '], for channel {} selected event'.format(k)
                h_1d = rt.TH1D(name, title, N_bins[0], -width+dx, width+dx)
                chain.Project(name, 'x_dut[{}]'.format(conf['idx_dut']), selection)
                h_1d.SetStats(0)
                h_1d.SetXTitle('x [mm]')
                h_1d.SetYTitle('Events / {:.1f} mm'.format(2*width/N_bins[0]))

                canvas['pos'][k] = rt.TCanvas('c_pos_'+str(k), 'c_pos_'+str(k), 800, 600)
                canvas['pos'][k].Divide(1,2)
                pad = canvas['pos'][k].cd(1)
                pad.SetRightMargin(0.15)
                h.DrawCopy('colz')
                pad = canvas['pos'][k].cd(2)
                pad.SetRightMargin(0.15)
                h_1d.DrawCopy()
                canvas['pos'][k].Update()
                canvas['pos'][k].SaveAs(out_dir + '/PositionXY_raw_ch{:02d}'.format(k)+figform)

            if (('PosSel' in configurations.plots) or ('PosSel+' in configurations.plots)) and not (conf['type'][:3] in ['tim', 'amp']):
                canvas['pos_sel'][k] = rt.TCanvas('c_pos_sel_'+str(k), 'c_pos_sel_'+str(k), 800, 600)
                canvas['pos_sel'][k].Divide(1,2)
                pad = canvas['pos_sel'][k].cd(1)
                pad.SetRightMargin(0.15)

                name = 'h_pos_sel'
                title = 'Events average selection efficiency'
                h = rt.TH3D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 11, -0.05, 1.05)

                chain.Project(name, selection + ':' + var, configurations.TracksCleaning['cuts'])
                h = h.Project3DProfile('yx')
                h.SetYTitle('y [mm]')
                h.SetXTitle('x [mm]')
                h.SetZTitle('Selection efficiency')
                h.DrawCopy('colz')

                name = 'h_pos_sel1D'
                title = 'Events average selection efficiency'
                h = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, 11, -0.05, 1.05)
                sel_trk = selection
                sel_now = ''
                if 'space_sel' in conf.keys():
                    sel_trk += ' && ' + conf['space_sel']['cut_y']
                    sel_now = conf['space_sel']['cut_y']
                var1D = var[var.find(':')+1:]
                chain.Project('h_pos_sel1D', selection + ':' + var1D, sel_now)
                h = h.ProfileX()
                h.SetYTitle('Selection efficiency')
                h.SetXTitle('x [mm]')
                h.SetStats(0)
                pad = canvas['pos_sel'][k].cd(2)
                pad.SetRightMargin(0.15)
                h.DrawCopy()

                canvas['pos_sel'][k].Update()
                canvas['pos_sel'][k].SaveAs(out_dir + '/PositionXY_sel_ch{:02d}'.format(k)+figform)



            if 'PosWeight' in configurations.plots and not (conf['type'][:3] in ['tim', 'amp']):
                h = rt.TH2D('h_amp_aux'+str(k), title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                chain.Project('h_amp_aux'+str(k), var, selection)

                ''' -------------Avg amp-----------'''
                name = 'h_amp_weight_pos_'+str(k)
                title = 'Average peak amplitude channel {} in selected event'.format(k)
                h_w = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h_w.SetXTitle('x [mm]')
                h_w.SetYTitle('y [mm]')
                h_w.SetZTitle('Average Amplitude [mV]')

                weights = '('+ selection +') * amp[' + str(k) + ']'
                chain.Project(name, var, weights)

                h_w.Divide(h)
                h_w.GetZaxis().SetRangeUser(h_w.GetMinimum(10), h_w.GetMaximum())
                h_w.SetStats(0)

                canvas['amp_w_pos'][k] = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 600)
                canvas['amp_w_pos'][k].Divide(1,2)
                pad = canvas['amp_w_pos'][k].cd(1)
                pad.SetRightMargin(0.15)
                h_w.DrawCopy('colz')

                pad = canvas['amp_w_pos'][k].cd(2)
                pad.SetRightMargin(0.15)
                h1D = rt.TProfile('h_aux', '', N_bins[0], -width_x+dx, width_x+dx)
                h1D.SetXTitle('x [mm]')
                h1D.SetYTitle('Amplitude [mV]')
                chain.Project('h_aux', 'amp[' + str(k) + ']' + ':' + var1D, sel_trk, 'prof')
                h1D.SetStats(0)
                h1D.DrawCopy()

                canvas['amp_w_pos'][k].Update()
                canvas['amp_w_pos'][k].SaveAs(out_dir + '/PositionXY_amp_weight_ch{:02d}'.format(k)+figform)

                ''' -------------Avg integral-----------'''
                name = 'h_int_weight_pos_'+str(k)
                title = 'Average integral channel {} in selected event'.format(k)
                h_w = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h_w.SetXTitle('x [mm]')
                h_w.SetYTitle('y [mm]')
                h_w.SetZTitle('Average integral [pC]')

                weights = '('+ selection +') * -integral[' + str(k) + ']'
                chain.Project(name, var, weights)

                h_w.Divide(h)
                arr, pos = rootTH2_to_np(h_w)
                arr = arr[arr!=0]
                h_w.GetZaxis().SetRangeUser(np.percentile(arr,10), np.percentile(arr,99))
                h_w.SetStats(0)

                canvas['int_w_pos'][k] = rt.TCanvas('c_int_w_pos_'+str(k), 'c_int_w_pos_'+str(k), 800, 600)
                canvas['int_w_pos'][k].Divide(1,2)
                pad = canvas['int_w_pos'][k].cd(1)
                pad.SetRightMargin(0.15)
                h_w.DrawCopy('colz')

                pad = canvas['int_w_pos'][k].cd(2)
                pad.SetRightMargin(0.15)
                h1D = rt.TProfile('h_aux', '', N_bins[0], -width_x+dx, width_x+dx)
                h1D.SetXTitle('x [mm]')
                h1D.SetYTitle('Integral [pC]')
                chain.Project('h_aux', '-integral[' + str(k) + ']' + ':' + var1D, sel_trk, 'prof')
                h1D.SetStats(0)
                h1D.DrawCopy()

                canvas['int_w_pos'][k].Update()
                canvas['int_w_pos'][k].SaveAs(out_dir + '/PositionXY_int_weight_ch{:02d}'.format(k)+figform)

                ''' -------------Avg risetime-----------'''
                name = 'h_risetime_weight_pos_'+str(k)
                title = 'Average risetime channel {} in selected event'.format(k)
                h_w = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h_w.SetXTitle('x [mm]')
                h_w.SetYTitle('y [mm]')
                h_w.SetZTitle('Average Slew rate [mV/ns]')

                weights = '('+ selection +') * -risetime[' + str(k) + ']'
                chain.Project(name, var, weights)

                h_w.Divide(h)
                arr, pos = rootTH2_to_np(h_w)
                arr = arr[arr!=0]
                h_w.GetZaxis().SetRangeUser(np.percentile(arr,10), np.percentile(arr,99))
                h_w.SetStats(0)
                canvas['int_w_pos'][k] = rt.TCanvas('c_risetime_w_pos_'+str(k), 'c_risetime_w_pos_'+str(k), 800, 600)
                canvas['int_w_pos'][k].Divide(1,2)
                pad = canvas['int_w_pos'][k].cd(1)
                pad.SetRightMargin(0.15)
                h_w.DrawCopy('colz')

                pad = canvas['int_w_pos'][k].cd(2)
                pad.SetRightMargin(0.15)
                h1D = rt.TProfile('h_aux', '', N_bins[0], -width_x+dx, width_x+dx)
                h1D.SetXTitle('x [mm]')
                h1D.SetYTitle('Slew rate [mV/ns]')
                chain.Project('h_aux', '-risetime[' + str(k) + ']' + ':' + var1D, sel_trk, 'prof')
                h1D.SetStats(0)
                h1D.DrawCopy()

                canvas['int_w_pos'][k].Update()
                canvas['int_w_pos'][k].SaveAs(out_dir + '/PositionXY_risetime_weight_ch{:02d}'.format(k)+figform)

        '''=========================== End Selections ==========================='''
        conf['sel'] =  selection
        if args.verbose:
            print 'selection:', conf['sel']
        # if 'space_sel' in conf.keys():
            # print conf['space_sel']['cut']

    '''End of single channels selection loop'''

    print '\n\n======================= Channel couples performances loop ==========================='

    headout_dir = out_dir
    best_result = {}
    for kL in configurations.ch_ordered:
        confL = configurations.channel[kL]
        #Look for left SiPM
        re_out = re.search(r'time([0-9]+)L', confL['type'])
        if (hasattr(re_out, 'group')) and confL['idx_ref'] >= 0:
            N_bar = int(re_out.group(1))
            found = False
            for kR in configurations.ch_ordered:
                confR = configurations.channel[kR]
                re_out = re.search(r'time([0-9]+)R', confR['type'])
                if (hasattr(re_out, 'group')):
                    if int(re_out.group(1)) == N_bar:
                        found =  True
                        break
            if not found:
                print 'Ch {} is left of bar {}, but no corresponding right found'.format(kL, N_bar)
                continue
        else:
            continue
        print '---> Bar #', N_bar, ': Left ch {} - Right ch {}'.format(kL, kR)

        results = []
        file_results = open(headout_dir+'/TimeResolution_bar{}.txt'.format(N_bar),'w')
        ln = '#avg_time_res, d_avg_time_res, var_time, var_ref, Correction\n'
        file_results.write(ln)

        best_result[N_bar] = Bauble()
        best_result[N_bar].dT = [-1, -1]
        best_result[N_bar].var = [None, None]
        best_result[N_bar].Correction = None

        if confL['idx_ref'] != confR['idx_ref']:
            print 'Inconsistency in reference chnnel'
        ref_type = configurations.channel[confL['idx_ref']]['type']

        BarInfo = Bauble()
        BarInfo.sel = [confL['sel'], confR['sel'], configurations.channel[confL['idx_ref']]['sel']]

        for kk in configurations.ch_ordered:
            aux_conf = configurations.channel[kk]['type']
            if 'amp'+str(kL) in aux_conf:
                confL['amp_ch'] = kk
                BarInfo.sel += [configurations.channel[kk]['sel']]
            if 'amp'+str(kR) in aux_conf:
                confR['amp_ch'] = kk
                BarInfo.sel += [configurations.channel[kk]['sel']]
            if not 'amp_ch' in confL.keys():
                    confL['amp_ch'] = kL
            if not 'amp_ch' in confR.keys():
                    confR['amp_ch'] = kR

        ''' ---------------  Position plots  ---------------'''
        var = 'y_dut[{}]:x_dut[{}]'.format(confL['idx_dut'], confL['idx_dut'])
        dy = configurations.xy_center[1]
        dx = configurations.xy_center[0]
        width_x = configurations.xy_center[2]
        width_y = configurations.xy_center[3]
        N_bins = [50, 50]

        if ('PosSel' in configurations.plots) or ('PosSel+' in configurations.plots):
            can = rt.TCanvas('c_pos_sel_'+str(N_bar), 'c_pos_sel_'+str(N_bar), 800, 500)
            can.Divide(2,1)
            can.cd(1)
            rt.gPad.Divide(1,3)
            pad = rt.gPad.cd(1)
            pad.SetRightMargin(0.15)

            name = 'h_pos_sel'
            title = 'Events average selection efficiency (R & L)'
            h = rt.TH3D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 11, -0.05, 1.05)

            chain.Project(name, ' && '.join(BarInfo.sel) + ':' + var, configurations.TracksCleaning['cuts'])
            h = h.Project3DProfile('yx')
            h.SetYTitle('y [mm]')
            h.SetXTitle('x [mm]')
            h.SetZTitle('Selection efficiency')
            h.SetStats(0)
            h.DrawCopy('colz')

            if ('PosSel+' in configurations.plots) and ('shape' in confL.keys() and 'shape' in confR.keys()):
                if confL['shape'] != confR['shape']:
                    print 'Shape right and left incongruent for bar', N_bar
                    continue

                if not confL['shape'] == 'None':
                    size = confL['shape'].split('x')

                    x = h.GetXaxis().GetBinCenter(1) - h.GetXaxis().GetBinWidth(1)*0.51 + float(size[0])
                    nbx = h.GetXaxis().FindBin(x)

                    y = h.GetYaxis().GetBinCenter(1) - h.GetYaxis().GetBinWidth(1)*0.51 + float(size[1])
                    nby = h.GetYaxis().FindBin(y)

                    arr_h, pos_h = rootTH2_to_np(h, cut=0.2)
                    # print arr_h

                    p_max = 0
                    idx_max = [0, 0]
                    for iy in range(0, arr_h.shape[0]-nby):
                        for ix in range(0, arr_h.shape[1]-nbx):
                            p = np.sum(arr_h[iy:iy+nby, ix:ix+nbx])
                            if p > p_max:
                                p_max = p
                                idx_max = [iy, ix]

                    avg_prob = p_max/(nbx*nby)
                    iy, ix = idx_max
                    x_start = h.GetXaxis().GetBinCenter(ix+1) - 0.5*h.GetXaxis().GetBinWidth(ix+1)
                    x_stop = h.GetXaxis().GetBinCenter(ix+nbx) + 0.5*h.GetXaxis().GetBinWidth(ix+nbx)
                    y_start = h.GetYaxis().GetBinCenter(iy+1) - 0.5*h.GetYaxis().GetBinWidth(iy+1)
                    y_stop = h.GetYaxis().GetBinCenter(iy+nby) + 0.5*h.GetYaxis().GetBinWidth(iy+nby)

                    best_theta = 0
                    cx = 0.5*(x_stop + x_start)
                    cy = 0.5*(y_stop + y_start)
                    x = '(x_dut[{}] - {})'.format(conf['idx_dut'], cx)
                    y = '(y_dut[{}] - {})'.format(conf['idx_dut'], cy)
                    var_template = '{c}*{y} - {s}*{x} + {cy}:{c}*{x} + {s}*{y} + {cx}'
                    for theta in np.arange(-0.1, 0.1, 0.005):
                        var_aux = var_template.format(x=x, y=y, c=np.cos(theta), s=np.sin(theta), cy=cy, cx=cx)

                        h = rt.TH3D('h_aux', title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 2, -0.5, 1.5)
                        chain.Project('h_aux', ' && '.join(BarInfo.sel) + ':' + var_aux, configurations.TracksCleaning['cuts'])
                        h = h.Project3DProfile('yx')
                        arr_h, pos_h = rootTH2_to_np(h, cut=0.2)
                        p = np.sum(arr_h[iy:iy+nby, ix:ix+nbx])
                        if p > p_max:
                            p_max = p
                            best_theta = theta

                    BarInfo.theta_z = best_theta
                    if args.verbose:
                        print 'Rotation of {:1.1e} rad detected aroud z'.format(best_theta)
                    BarInfo.x_rot = '{c}*{x} + {s}*{y} + {cx}'.format(x=x, y=y, c=np.cos(best_theta), s=np.sin(best_theta), cx=cx)
                    BarInfo.y_rot = '{c}*{y} - {s}*{x} + {cy}'.format(x=x, y=y, c=np.cos(best_theta), s=np.sin(best_theta), cy=cy)
                    if args.verbose:
                        print 'x_var:', BarInfo.x_rot
                        print 'y_var:', BarInfo.y_rot

                    var = BarInfo.y_rot+':'+BarInfo.x_rot
                    h = rt.TH3D('h_aux', title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 2, -0.5, 1.5)
                    chain.Project('h_aux', ' && '.join(BarInfo.sel) + ':' + var, configurations.TracksCleaning['cuts'])
                    h = h.Project3DProfile('yx')
                    h.SetYTitle('y [mm]')
                    h.SetXTitle('x [mm]')
                    h.SetZTitle('Selection efficiency')
                    h.SetStats(0)
                    h.DrawCopy('colz')

                    line.SetLineStyle(9)
                    line.SetLineColor(6)
                    line.SetLineWidth(2)
                    line.DrawLine(x_start, y_start, x_stop, y_start)
                    line.DrawLine(x_start, y_stop, x_stop, y_stop)
                    line.DrawLine(x_start, y_start, x_start, y_stop)
                    line.DrawLine(x_stop, y_start, x_stop, y_stop)

                    BarInfo.space_limits = [x_start, x_stop, y_start, y_stop]
                    BarInfo.cut_y = '({y} < {h} && {y} > {l})'.format(y=BarInfo.y_rot, l=y_start, h=y_stop)
                    BarInfo.cut_x = '({x} < {h} && {x} > {l})'.format(x=BarInfo.x_rot, l=x_start, h=x_stop)
                    if args.verbose:
                        print 'Space limits:', BarInfo.space_limits

            can.cd(1)
            pad = rt.gPad.cd(2)
            pad.SetRightMargin(0.15)
            name = 'h_pos_selL'
            title = 'Events average selection efficiency (Left only - ch {})'.format(kL)
            h = rt.TH3D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 11, -0.05, 1.05)
            chain.Project(name, confL['sel'] + ':' + var, configurations.TracksCleaning['cuts'])
            h = h.Project3DProfile('yx')
            h.SetYTitle('y [mm]')
            h.SetXTitle('x [mm]')
            h.SetZTitle('Selection efficiency')
            h.SetStats(0)
            h.DrawCopy('colz')

            can.cd(1)
            pad = rt.gPad.cd(3)
            pad.SetRightMargin(0.15)
            name = 'h_pos_selR'
            title = 'Events average selection efficiency (Right only - ch {})'.format(kR)
            h = rt.TH3D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy, 11, -0.05, 1.05)
            chain.Project(name, confR['sel'] + ':' + var, configurations.TracksCleaning['cuts'])
            h = h.Project3DProfile('yx')
            h.SetYTitle('y [mm]')
            h.SetXTitle('x [mm]')
            h.SetZTitle('Selection efficiency')
            h.SetStats(0)
            h.DrawCopy('colz')

            pad = can.cd(2)
            leg = rt.TLegend(0.85,0.8,0.98,0.93)

            name = 'h_pos_sel1D'
            title = 'Events average selection efficiency within bar area'
            h = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, 11, -0.05, 1.05)
            sel_trk = ' && '.join(BarInfo.sel)
            sel_now = configurations.TracksCleaning['cuts']
            if hasattr(BarInfo, 'cut_y'):
                sel_now += ' && ' + BarInfo.cut_y
            var1D = var[var.find(':')+1:]
            chain.Project('h_pos_sel1D', sel_trk + ':' + var1D, sel_now)
            h = h.ProfileX()
            h.SetYTitle('Selection efficiency')
            h.SetXTitle('x [mm]')
            h.SetStats(0)
            h.SetLineColor(1)
            h.GetYaxis().SetRangeUser(0., 1.)
            leg.AddEntry(h, 'L+R', 'l')
            h.DrawCopy()

            h = rt.TH2D(name+'L', title, N_bins[0], -width_x+dx, width_x+dx, 11, -0.05, 1.05)
            sel_trk = [BarInfo.sel[0], BarInfo.sel[2]]
            if len(BarInfo.sel) > 3:
                sel_trk += [configurations.channel[confL['amp_ch']]['sel']]
            sel_trk = ' && '.join(sel_trk)
            var1D = var[var.find(':')+1:]
            chain.Project('h_pos_sel1DL', sel_trk + ':' + var1D, sel_now)
            hL = h.ProfileX()
            hL.SetStats(0)
            hL.SetLineColor(2)
            leg.AddEntry(hL, 'L', 'l')
            hL.DrawCopy('SAME')

            h = rt.TH2D(name+'R', title, N_bins[0], -width_x+dx, width_x+dx, 11, -0.05, 1.05)
            sel_trk = [BarInfo.sel[1], BarInfo.sel[2]]
            if len(BarInfo.sel) > 3:
                sel_trk += [configurations.channel[confR['amp_ch']]['sel']]
            sel_trk = ' && '.join(sel_trk)
            var1D = var[var.find(':')+1:]
            chain.Project('h_pos_sel1DR', sel_trk + ':' + var1D, sel_now)
            hR = h.ProfileX()
            hR.SetStats(0)
            hR.SetLineColor(4)
            leg.AddEntry(hR, 'R', 'l')
            hR.DrawCopy('SAME')

            leg.Draw('SAME')

            can.Update()
            can.SaveAs(out_dir + '/PositionXY_sel_Bar{:02d}'.format(N_bar)+figform)

        if 'PosWeight' in configurations.plots:
            h = rt.TH2D('h_amp_aux'+str(k), title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            chain.Project('h_amp_aux'+str(k), var, ' && '.join(BarInfo.sel))

            ''' -------------Avg amp-----------'''
            name = 'h_amp_weight_pos_'+str(N_bar)
            title = 'Average peak amplitude (L+R) bar {} in selected event'.format(N_bar)
            h2D_sum = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_sum.SetXTitle('x [mm]')
            h2D_sum.SetYTitle('y [mm]')
            h2D_sum.SetZTitle('Average Amplitude [mV]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * (amp[{}] + amp[{}])'.format(confL['amp_ch'], confR['amp_ch'])
            chain.Project(name, var, weights)
            h2D_sum.Divide(h)
            h2D_sum.GetZaxis().SetRangeUser(h2D_sum.GetMinimum(10), h2D_sum.GetMaximum())
            h2D_sum.SetStats(0)

            name = 'h_ampL_weight_pos_'+str(N_bar)
            title = 'Average peak amplitude (L only - ch {}) in selected event'.format(confL['amp_ch'])
            h2D_L = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_L.SetXTitle('x [mm]')
            h2D_L.SetYTitle('y [mm]')
            h2D_L.SetZTitle('Average Amplitude [mV]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * amp[{}]'.format(confL['amp_ch'])
            chain.Project(name, var, weights)
            h2D_L.Divide(h)
            h2D_L.GetZaxis().SetRangeUser(h2D_L.GetMinimum(10), h2D_L.GetMaximum())
            h2D_L.SetStats(0)

            name = 'h_ampR_weight_pos_'+str(N_bar)
            title = 'Average peak amplitude (R only - ch {}) in selected event'.format(confR['amp_ch'])
            h2D_R = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_R.SetXTitle('x [mm]')
            h2D_R.SetYTitle('y [mm]')
            h2D_R.SetZTitle('Average Amplitude [mV]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * amp[{}]'.format(confR['amp_ch'])
            chain.Project(name, var, weights)
            h2D_R.Divide(h)
            h2D_R.GetZaxis().SetRangeUser(h2D_R.GetMinimum(10), h2D_R.GetMaximum())
            h2D_R.SetStats(0)


            sel_1D = ' && '.join(BarInfo.sel)
            if hasattr(BarInfo, 'cut_y'):
                sel_1D += ' && ' + BarInfo.cut_y

            h1D_sum = rt.TProfile('h_sum', 'Average peak amplitude in selected events', N_bins[0], -width_x+dx, width_x+dx)
            h1D_sum.SetXTitle('x [mm]')
            h1D_sum.SetYTitle('Amplitude [mV]')
            var1D = var[var.find(':')+1:]
            chain.Project('h_sum', 'amp[{}] + amp[{}]'.format(confL['amp_ch'], confR['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_sum.SetStats(0)

            h1D_L = rt.TProfile('h_L', 'Average peak amplitude in selected events', N_bins[0], -width_x+dx, width_x+dx)
            chain.Project('h_L', 'amp[{}]'.format(confL['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_L.SetStats(0)

            h1D_R = rt.TProfile('h_R', 'Average peak amplitude in selected events', N_bins[0], -width_x+dx, width_x+dx)
            chain.Project('h_R', 'amp[{}]'.format(confR['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_R.SetStats(0)

            can = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 500)
            can.Divide(2,1)
            padL = can.cd(1)
            padL.Divide(1,3)

            pad = padL.cd(1)
            pad.SetRightMargin(0.15)
            h2D_sum.DrawCopy('colz')
            if hasattr(BarInfo, 'space_limits'):
                x_start, x_stop, y_start, y_stop = BarInfo.space_limits
                line.SetLineStyle(9)
                line.SetLineColor(6)
                line.SetLineWidth(2)
                line.DrawLine(x_start, y_start, x_stop, y_start)
                line.DrawLine(x_start, y_stop, x_stop, y_stop)
                line.DrawLine(x_start, y_start, x_start, y_stop)
                line.DrawLine(x_stop, y_start, x_stop, y_stop)

            pad = padL.cd(2)
            pad.SetRightMargin(0.15)
            h2D_L.DrawCopy('colz')
            pad = padL.cd(3)
            pad.SetRightMargin(0.15)
            h2D_R.DrawCopy('colz')


            pad = can.cd(2)
            leg = rt.TLegend(0.85,0.8,0.95,0.92)

            h1D_sum.SetLineColor(1)
            h1D_sum.GetYaxis().SetRangeUser(0., 1.2*h1D_sum.GetMaximum())
            h1D_sum.Draw()
            leg.AddEntry(h1D_sum, 'L+R', 'l')
            h1D_L.SetLineColor(2)
            h1D_L.Draw('SAME')
            leg.AddEntry(h1D_L, 'L', 'l')
            h1D_R.SetLineColor(4)
            h1D_R.Draw('SAME')
            leg.AddEntry(h1D_R, 'R', 'l')
            leg.Draw()

            can.Update()
            can.SaveAs(out_dir + '/PositionXY_amp_weight_bar{:02d}'.format(N_bar)+figform)

            AmpAvg_sum= np.mean([h1D_sum.GetBinContent(i)  for i in range(h1D_sum.GetNbinsX()) if (h1D_sum.GetBinCenter(i)>x_start and h1D_sum.GetBinCenter(i)<x_stop)])
            AmpAvg_sum_err= np.mean([h1D_sum.GetBinError(i)  for i in range(h1D_sum.GetNbinsX()) if (h1D_sum.GetBinCenter(i)>x_start and h1D_sum.GetBinCenter(i)<x_stop)])
            AmpAvg_L= np.mean([h1D_L.GetBinContent(i) for i in range(h1D_sum.GetNbinsX()) if (h1D_L.GetBinCenter(i)>x_start and h1D_L.GetBinCenter(i)<x_stop) ])
            AmpAvg_L_err= np.mean([h1D_L.GetBinError(i)  for i in range(h1D_L.GetNbinsX()) if (h1D_L.GetBinCenter(i)>x_start and h1D_L.GetBinCenter(i)<x_stop)])
            AmpAvg_R= np.mean([h1D_R.GetBinContent(i) for i in range(h1D_R.GetNbinsX()) if (h1D_R.GetBinCenter(i)>x_start and h1D_R.GetBinCenter(i)<x_stop) ])
            AmpAvg_R_err= np.mean([h1D_R.GetBinError(i)  for i in range(h1D_R.GetNbinsX()) if (h1D_R.GetBinCenter(i)>x_start and h1D_R.GetBinCenter(i)<x_stop)])

            file_ampres = open(headout_dir+'/Amplitude_bar{}.txt'.format(N_bar),'w')
            ln = '#avg_amp_LplusR, avg_amp_LplusR_err, avg_amp_L, avg_amp_L_err, avg_amp_R, avg_amp_R_err\n'
            file_ampres.write(ln)
            ln = '{:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}\n'.format(AmpAvg_sum,AmpAvg_sum_err,AmpAvg_L,AmpAvg_L_err,AmpAvg_R,AmpAvg_R_err)
            file_ampres.write(ln)
            file_ampres.close()
            print "Average amplitude : ", AmpAvg_sum, "+/- ", AmpAvg_sum_err
            print "Average left amplitude : ", AmpAvg_L, "+/- ", AmpAvg_L_err
            print "Average right amplitude : ", AmpAvg_R, "+/- ", AmpAvg_R_err

            ''' -------------Avg integral-----------'''

            name = 'h_weight_pos_'+str(N_bar)
            title = 'Average integral (L+R) bar {} in selected event'.format(N_bar)
            h2D_sum = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_sum.SetXTitle('x [mm]')
            h2D_sum.SetYTitle('y [mm]')
            h2D_sum.SetZTitle('Average integral [pC]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * (-integral[{}] -integral[{}])'.format(confL['amp_ch'], confR['amp_ch'])
            chain.Project(name, var, weights)
            h2D_sum.Divide(h)
            arr, pos = rootTH2_to_np(h2D_sum)
            arr = arr[arr!=0]
            h2D_sum.GetZaxis().SetRangeUser(h2D_sum.GetMinimum(10), h2D_sum.GetMaximum())
            h2D_sum.SetStats(0)

            name = 'h_L_weight_pos_'+str(N_bar)
            title = 'Average integral (L only - ch {}) in selected event'.format(confL['amp_ch'])
            h2D_L = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_L.SetXTitle('x [mm]')
            h2D_L.SetYTitle('y [mm]')
            h2D_L.SetZTitle('Integral [pC]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * -integral[{}]'.format(confL['amp_ch'])
            chain.Project(name, var, weights)
            h2D_L.Divide(h)
            arr, pos = rootTH2_to_np(h2D_L)
            arr = arr[arr!=0]
            h2D_L.GetZaxis().SetRangeUser(h2D_L.GetMinimum(10), h2D_L.GetMaximum())
            h2D_L.SetStats(0)

            name = 'h_R_weight_pos_'+str(N_bar)
            title = 'Average integral (R only - ch {}) in selected event'.format(confR['amp_ch'])
            h2D_R = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
            h2D_R.SetXTitle('x [mm]')
            h2D_R.SetYTitle('y [mm]')
            h2D_R.SetZTitle('Integral [pC]')
            weights = '('+ ' && '.join(BarInfo.sel) +') * -integral[{}]'.format(confR['amp_ch'])
            chain.Project(name, var, weights)
            h2D_R.Divide(h)
            arr, pos = rootTH2_to_np(h2D_R)
            arr = arr[arr!=0]
            h2D_R.GetZaxis().SetRangeUser(h2D_R.GetMinimum(10), h2D_R.GetMaximum())
            h2D_R.SetStats(0)


            h1D_sum = rt.TProfile('h_sum', 'Average integral in selected events', N_bins[0], -width_x+dx, width_x+dx)
            h1D_sum.SetXTitle('x [mm]')
            h1D_sum.SetYTitle('Integral [pC]')
            chain.Project('h_sum', '-integral[{}] -integral[{}]'.format(confL['amp_ch'], confR['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_sum.SetStats(0)

            h1D_L = rt.TProfile('h_L', 'Average integral in selected events', N_bins[0], -width_x+dx, width_x+dx)
            chain.Project('h_L', '-integral[{}]'.format(confL['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_L.SetStats(0)

            h1D_R = rt.TProfile('h_R', 'Average integral in selected events', N_bins[0], -width_x+dx, width_x+dx)
            chain.Project('h_R', '-integral[{}]'.format(confR['amp_ch']) + ':' + var1D, sel_1D, 'prof')
            h1D_R.SetStats(0)

            can = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 500)
            can.Divide(2,1)
            padL = can.cd(1)
            padL.Divide(1,3)

            pad = padL.cd(1)
            pad.SetRightMargin(0.15)
            h2D_sum.DrawCopy('colz')
            if hasattr(BarInfo, 'space_limits'):
                x_start, x_stop, y_start, y_stop = BarInfo.space_limits
                line.SetLineStyle(9)
                line.SetLineColor(6)
                line.SetLineWidth(2)
                line.DrawLine(x_start, y_start, x_stop, y_start)
                line.DrawLine(x_start, y_stop, x_stop, y_stop)
                line.DrawLine(x_start, y_start, x_start, y_stop)
                line.DrawLine(x_stop, y_start, x_stop, y_stop)

            pad = padL.cd(2)
            pad.SetRightMargin(0.15)
            h2D_L.DrawCopy('colz')
            pad = padL.cd(3)
            pad.SetRightMargin(0.15)
            h2D_R.DrawCopy('colz')


            pad = can.cd(2)
            leg = rt.TLegend(0.85,0.8,0.95,0.92)

            h1D_sum.SetLineColor(1)
            h1D_sum.GetYaxis().SetRangeUser(0., 1.2*h1D_sum.GetMaximum())
            h1D_sum.Draw()
            leg.AddEntry(h1D_sum, 'L+R', 'l')
            h1D_L.SetLineColor(2)
            h1D_L.Draw('SAME')
            leg.AddEntry(h1D_L, 'L', 'l')
            h1D_R.SetLineColor(4)
            h1D_R.Draw('SAME')
            leg.AddEntry(h1D_R, 'R', 'l')
            leg.Draw()

            can.Update()
            can.SaveAs(out_dir + '/PositionXY_int_weight_bar{:02d}'.format(N_bar)+figform)

        if hasattr(BarInfo, 'space_limits'):
            BarInfo.sel += [BarInfo.cut_y, BarInfo.cut_x]

        var_y = var[:var.find(':')]
        var_x = var[var.find(':')+1:]

        for v_ref, v_time in itertools.product(configurations.Var[ref_type], configurations.Var['time']):
            print 'Running on:', v_time, v_ref

            time_var_chref = v_ref+'[{}]'.format(confL['idx_ref'])
            time_var = '0.5*({v}[{L}]+{v}[{R}])'.format(v=v_time, L=kL, R=kR)
            spacelike_var = '({v}[{L}]-{v}[{R}])'.format(v=v_time, L=kL, R=kR)

            out_dir = '{}/{}_{}'.format(headout_dir, v_time, v_ref)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                shutil.copy('lib/index.php', out_dir + '/index.php')

            var_dT = time_var + ' - ' + time_var_chref
            var_dTL = '{v}[{L}] - '.format(v=v_time, L=kL) + time_var_chref
            var_dTR = '{v}[{R}] - '.format(v=v_time, R=kR) + time_var_chref

            if args.verbose:
                print 'var_dT = \'{}\''.format(var_dT)
                print 'var_dTL = \'{}\''.format(var_dTL)
                print 'var_dTR = \'{}\''.format(var_dTR)

            selection = BarInfo.sel
            selection += ['({v}[{L}]!= 0 && {v}[{R}]!=0)'.format(v=v_time, L=kL, R=kR)]
            selection = ' && '.join(selection)
            print selection
            if chain.GetEntries(selection) < 5:
                print chain.GetEntries(selection)
                print 'Not enought stat ({})'.format(chain.GetEntries(selection))
                continue


            '''=========================== Space like variable ==========================='''
            if 'SpaceLike' in configurations.plots:
                aux_x =  tree2array(chain, BarInfo.x_rot, selection).flatten()
                aux_x -= np.min(aux_x)
                aux_spacelike =  tree2array(chain, spacelike_var, selection).flatten()

                inputs = np.column_stack((np.ones_like(aux_x), aux_x))
                coeff, r, rank, s = np.linalg.lstsq(inputs, aux_spacelike)
                # coeff, r, rank, s = np.linalg.lstsq(inputs, aux_spacelike, rcond=None)
                b = [None, None, None, None, np.percentile(aux_spacelike, 1), np.percentile(aux_spacelike, 99)]
                h_SLvsX = create_TH2D(np.column_stack((aux_x, aux_spacelike)),
                                     'h_SLvsX', 'Bar '+str(N_bar),
                                     binning=b,
                                     axis_title=['x [mm]', spacelike_var + ' [ns]']
                                     )
                h_SLvsX.SetStats(0)
                canvas['dt_vs_amp'][k] = rt.TCanvas('h_SLvsX'+str(N_bar), 'h_SLvsX'+str(N_bar), 800, 600)
                h_SLvsX.DrawCopy('colz')
                prof = h_SLvsX.ProfileX('prof_amp')
                prof.SetLineColor(6)
                prof.SetLineWidth(2)
                prof.DrawCopy('SAMEE1')

                f = rt.TF1('SLvsX_fit'+str(k),'[0]+[1]*x', np.min(aux_x), np.max(aux_x))
                for j,a in enumerate(coeff):
                    f.SetParameter(j, a)
                f.SetLineColor(6)
                f.SetLineStyle(9)
                f.DrawCopy('SAMEL')

                l = rt.TLatex()
                l.SetTextSize(0.04)
                v = 2./coeff[1]
                length = -coeff[0]*v
                v2 = -50./coeff[0]
                # v3 = np.sum(np.square(2*aux_x - 50.))/np.sum((2*aux_x - 50.)*aux_spacelike)
                # print 'v3 = {:.0f} mm/ns'.format(v3)
                msg = 'v = {:.0f} ({:.0f}) mm/ns - l = {:.1f} mm'.format(v, v2, length)
                print msg
                print 'Refractive index = {:.2f}'.format(299.0/v2)
                l.DrawLatexNDC(0.15, 0.89, msg)

                canvas['dt_vs_amp'][k].Update()
                canvas['dt_vs_amp'][k].SaveAs(out_dir + '/SLvsX_{}'.format(N_bar)+figform)


            '''=========================== Raw time resolution ==========================='''
            if 'TimeResRaw' in configurations.plots:
                delta_t =  tree2array(chain, var_dT, selection).flatten()
                if ( len(delta_t) ==0):
                    print 'Empty delta'
                    continue

                if os.uname()[1].startswith('lxplus'):
                    median = np.percentile(delta_t, 50)[0]
                    width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))[0]
                else:
                    median = np.percentile(delta_t, 50)
                    width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))

                name = 'h_delta_t_raw_'+str(k)
                if width == 0:
                    width = np.std(delta_t)
                if width == 0 or np.isnan(width):
                    width = 0.1
                if np.isnan(median):
                    continue
                title = 'Time resolution for channel '+str(k)+', width 10-90 = {:.2f} ns'.format(width)
                h = create_TH1D(delta_t, name, title,
                                    binning = [ None, median-2*width, median+2*width],
                                    axis_title = [var_dT + ' [ns]', 'Events'])

                rt.gStyle.SetStatX(1.)
                canvas['t_res_raw'][N_bar] = rt.TCanvas('c_t_res_raw_'+str(k), 'c_t_res_raw_'+str(k), 800, 600)
                h.Fit('gaus', 'LQR','', median-width, median+width)
                h = h.DrawCopy('LE')
                canvas['t_res_raw'][N_bar].Update()
                canvas['t_res_raw'][N_bar].SaveAs(out_dir + '/TimeResolution_raw_bar{:02d}'.format(N_bar)+figform)

                selection += '&& ({v} < {h} && {v} > {l})'.format(v=var_dT, h=median+2*width, l=median-2*width)
                if args.verbose:
                    print 'selection:', selection

            '''=========================== Mean time of arrival ==========================='''
            if 'MeanTime' in configurations.plots:
                h = rt.TH2D('h_amp_aux'+str(k), title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                chain.Project('h_amp_aux'+str(k), var, selection)

                ''' -------------Time of arrival-----------'''
                name = 'h_aux_'+str(N_bar)
                title = 'Average time stamp (L+R)/2 - bar {} in selected event'.format(N_bar)
                h2D_sum = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h2D_sum.SetXTitle('x [mm]')
                h2D_sum.SetYTitle('y [mm]')
                h2D_sum.SetZTitle('Average time [ns]')
                weights = '('+ selection +') * ('+var_dT+')'
                chain.Project(name, var, weights)
                h2D_sum.Divide(h)
                h2D_sum.GetZaxis().SetRangeUser(median-2*width, median+2*width)
                h2D_sum.SetStats(0)

                name = 'h_auxL_'+str(N_bar)
                title = 'Average mean time stamp Left - bar {} in selected event'.format(N_bar)
                h2D_L = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h2D_L.SetXTitle('x [mm]')
                h2D_L.SetYTitle('y [mm]')
                h2D_L.SetZTitle('Average time [ns]')
                weights = '('+ selection +') * ('+var_dTL+')'
                chain.Project(name, var, weights)
                h2D_L.Divide(h)
                h2D_L.GetZaxis().SetRangeUser(median-2*width, median+2*width)
                h2D_L.SetStats(0)

                name = 'h_auxR_'+str(N_bar)
                title = 'Average time stamp Right - bar {} in selected event'.format(N_bar)
                h2D_R = rt.TH2D(name, title, N_bins[0], -width_x+dx, width_x+dx, N_bins[1], -width_y+dy, width_y+dy)
                h2D_R.SetXTitle('x [mm]')
                h2D_R.SetYTitle('y [mm]')
                h2D_R.SetZTitle('Average time [ns]')
                weights = '('+ selection +') * ('+var_dTR+')'
                chain.Project(name, var, weights)
                h2D_R.Divide(h)
                h2D_R.GetZaxis().SetRangeUser(median-2*width, median+2*width)
                h2D_R.SetStats(0)

                h1D_sum = rt.TProfile('h_sum', 'Average time stamp in selected events', N_bins[0], -width_x+dx, width_x+dx)
                h1D_sum.SetXTitle('x [mm]')
                h1D_sum.SetYTitle('Average time [ns]')
                chain.Project('h_sum', var_dT + ':' + var_x, selection, 'prof')
                h1D_sum.SetStats(0)
                h1D_sum.GetYaxis().SetRangeUser(median-2*width, median+2*width)

                h1D_L = rt.TProfile('h_L', '', N_bins[0], -width_x+dx, width_x+dx)
                chain.Project('h_L', var_dTL + ':' + var_x, selection, 'prof')
                h1D_L.SetStats(0)

                h1D_R = rt.TProfile('h_R', '', N_bins[0], -width_x+dx, width_x+dx)
                chain.Project('h_R', var_dTR + ':' + var_x, selection, 'prof')
                h1D_R.SetStats(0)

                Print2D_1D_canvas([h2D_sum, h2D_L, h2D_R], [h1D_sum, h1D_L, h1D_R], out_dir + '/TimeMean_raw_bar{:02d}'.format(N_bar)+figform)

            '''=========================== Time resolution 2D ==========================='''
            if 'TimeResRaw2D' in configurations.plots:
                if not (('shape' in confL.keys()) and ('shape' in confR.keys())):
                    print 'Please specify a shape for both ends in order to print TimeResRaw2D'
                    continue

                b = [var_dT, BarInfo.x_rot, BarInfo.y_rot]
                if 'TimeRes2DAmp' in configurations.plots:
                    b += ['amp[{}]'.format(confL['amp_ch']), 'amp[{}]'.format(confR['amp_ch'])]
                b += [spacelike_var]

                data_raw = tree2array(chain, b, selection).view(np.recarray)
                data = np.zeros((data_raw.shape[0],3))
                if 'TimeRes2DAmp' in configurations.plots:
                    data = np.zeros((data_raw.shape[0],6))
                    for iev, d in enumerate(data_raw):
                        data[iev] = [d[0][0], d[1], d[2], d[3], d[4], d[5]]
                else:
                    for iev, d in enumerate(data_raw):
                        data[iev] = [d[0][0], d[1], d[2]]

                if ( len(data) ==0):
                    print 'Empty delta'
                    continue

                bin_width = configurations.space_bin_size

                lim = np.array(BarInfo.space_limits).astype(np.float)
                size = np.array(confL['shape'].split('x')).astype(np.float)

                x_sec, x_step = np.linspace(lim[0], lim[1], 1+int(size[0]/bin_width), retstep=True)
                y_sec, y_step = np.linspace(lim[2], lim[3], 1+int(size[1]/bin_width), retstep=True)

                if lim[1]-lim[0]+2*x_step >= lim[3]-lim[2]+2*y_step:
                    b = [x_sec.shape[0]+1, lim[0]-x_step, lim[1]+x_step]
                    space_diff = lim[1]-lim[0]+2*x_step - (lim[3]-lim[2])
                    adding_steps = int(space_diff/y_step)
                    if adding_steps%2 == 1:
                        adding_steps +=1
                    b += [y_sec.shape[0]-1+adding_steps, lim[2]-y_step*adding_steps/2, lim[3]+y_step*adding_steps/2]
                else:
                    b = [y_sec.shape[0]+1, lim[2]-y_step, lim[3]+y_step]
                    space_diff = lim[3]-lim[2]+2*y_step - (lim[1]-lim[0])
                    adding_steps = int(space_diff/x_step)
                    if adding_steps%2 == 1:
                        adding_steps +=1
                    b = [x_sec.shape[0]-1+adding_steps, lim[0]-x_step*adding_steps/2, lim[1]+x_step*adding_steps/2] + b

                name = 'h_TimeResRaw2D_bar{:02d}'.format(N_bar)
                title = 'Raw time res breakdown ch{:02d}'.format(k)
                h_2D_res_raw = rt.TH2D(name, title, b[0], b[1], b[2], b[3], b[4], b[5])
                ResRaw = []

                if 'TimeRes2DAmp' in configurations.plots:
                    name = 'h_TimeRes2DAmp_bar{:02d}'.format(N_bar)
                    title = 'Time resolution amplitude corrected breakdown bar{:02d}'.format(N_bar)
                    h_2D_res_amp = rt.TH2D(name, title, b[0], b[1], b[2], b[3], b[4], b[5])
                    ResAmp = []


                it = itertools.product(zip(x_sec[:-1], x_sec[1:]), zip(y_sec[:-1], y_sec[1:]))
                if not os.path.exists(out_dir + '/TimeResRaw2D_bar{:02d}'.format(N_bar)):
                    os.mkdir(out_dir + '/TimeResRaw2D_bar{:02d}'.format(N_bar))
                    shutil.copy('lib/index.php', out_dir + '/TimeResRaw2D_bar{:02d}'.format(N_bar)+'/index.php')
                if 'TimeRes2DAmp' in configurations.plots and not os.path.exists(out_dir + '/TimeRes2DAmp_bar{:02d}'.format(N_bar)):
                    os.mkdir(out_dir + '/TimeRes2DAmp_bar{:02d}'.format(N_bar))
                    shutil.copy('lib/index.php', out_dir +  '/TimeRes2DAmp_bar{:02d}'.format(N_bar)+'/index.php')

                for ib, ((xd, xu), (yd, yu)) in enumerate(it):
                    selx = np.logical_and(data[:,1]>xd, data[:,1]<xu)
                    sely = np.logical_and(data[:,2]>yd, data[:,2]<yu)
                    sel = np.logical_and(selx, sely)
                    if np.sum(sel) < 10:
                        continue
                    aux_d = data[sel]
                    dt = aux_d[:,0]

                    fill_TimeResHisto(N_bar, dt,
                                      h_2D_res_raw,
                                      ResRaw, 'TimeResRaw2D',
                                      'Raw',
                                      canvas,
                                      not ('TimeRes2DAmp' in configurations.plots))


                    if 'TimeRes2DAmp' in configurations.plots:
                        dt_sel = np.logical_and(dt>median-2*width, dt<median+2*width)
                        ampL = aux_d[:,3]
                        ampR = aux_d[:,4]
                        x_bar = aux_d[:,5]
                        amp_ratio = ampL/ampR
                        # inputs = np.column_stack((np.ones_like(amp_ratio), amp_ratio, amp_ratio**2))
                        inputs = np.column_stack((np.ones_like(ampL), ampL, ampR, np.square(ampL), ampR*ampL, np.square(ampR)))
                        coeff, r, rank, s = np.linalg.lstsq(inputs[dt_sel], dt[dt_sel])
                        # coeff, r, rank, s = np.linalg.lstsq(inputs[dt_sel], dt[dt_sel], rcond=None)
                        b = [None, None, None, None, median-2*width, median+2*width]
                        h_TvsA = create_TH2D(np.column_stack((amp_ratio, dt)),
                                             'h_TvsA', 'Bar '+str(N_bar),
                                             binning=b,
                                             axis_title=['Amp{} (L)/ Amp{} (R)'.format(kL, kR), var_dT + ' [ns]', '']
                                             )
                        canvas['dt_vs_amp'][k] = rt.TCanvas('dt_vs_amp'+str(k), 'dt_vs_amp'+str(k), 800, 600)
                        h_TvsA.DrawCopy('colz')
                        prof = h_TvsA.ProfileX('prof_amp')
                        prof.SetLineColor(6)
                        prof.SetLineWidth(2)
                        prof.DrawCopy('SAMEE1')

                        # f = rt.TF1('amp_fit'+str(k),'[0]+[1]*x+[2]*x^2', np.min(amp), np.max(amp))
                        # for j,a in enumerate(coeff):
                        #     f.SetParameter(j, a)
                        # f.SetLineColor(6)
                        # f.SetLineStyle(9)
                        # f.DrawCopy('SAMEL')
                        canvas['dt_vs_amp'][k].Update()
                        canvas['dt_vs_amp'][k].SaveAs(out_dir +  '/TimeRes2DAmp_bar{:02d}'.format(N_bar)+'/TvsAmp_{}'.format(ib)+figform)

                        dt_corr = dt - np.dot(inputs, coeff)

                        fill_TimeResHisto(N_bar, dt_corr,
                                          h_2D_res_amp,
                                          ResAmp, 'TimeRes2DAmp',
                                          'Amplitude corrected',
                                          canvas)


                ResRaw = np.array(ResRaw)
                tag = 'TimeResRaw2D_bar{:02d}'.format(N_bar)
                canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
                Set_2D_colz_graphics()
                h_2D_res_raw.GetZaxis().SetRangeUser(0.9*np.min(ResRaw[:,0]), 1.1*np.max(ResRaw[:,0]))
                h_2D_res_raw.SetStats(0)
                h_2D_res_raw.SetYTitle('y [mm]')
                h_2D_res_raw.SetXTitle('x [mm]')
                h_2D_res_raw.SetZTitle(var_dT+' [ps]')
                h_2D_res_raw.DrawCopy('colz')
                rt.gStyle.SetPaintTextFormat(".0f");
                h_2D_res_raw.DrawCopy('TEXT SAME')

                avg_time_res = np.average(ResRaw[:,0], weights=1./np.square(ResRaw[:,1]))
                d_avg_time_res = 1./np.sqrt(np.sum(1./np.square(ResRaw[:,1])))
                results.append([avg_time_res, d_avg_time_res, v_time, v_ref, False])
                ln = '{:.2f}  {:.2f}  {}  {}  {}\n'.format(avg_time_res, d_avg_time_res, v_time, v_ref, False)
                file_results.write(ln)
                if best_result[N_bar].dT[0] < 0 or  best_result[N_bar].dT[0] > avg_time_res:
                    best_result[N_bar].dT = [avg_time_res, d_avg_time_res]
                    best_result[N_bar].var = [v_time, v_ref]
                    best_result[N_bar].AmpCorr = False

                msg = 'Avg raw resolution: {:.1f} +/- {:.1f} ps'.format(avg_time_res, d_avg_time_res)
                print msg
                l = rt.TLatex()
                l.SetTextSize(0.04);
                l.DrawLatexNDC(0.15, 0.89, msg.replace('+/-', '#pm'))
                canvas['c_'+tag].Update()
                canvas['c_'+tag].SaveAs(out_dir + '/TimeResRaw2D_bar{:02d}'.format(N_bar)+figform)

                if 'TimeRes2DAmp' in configurations.plots:
                    ResAmp = np.array(ResAmp)
                    tag = 'TimeRes2DAmp_bar{:02d}'.format(N_bar)
                    canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
                    Set_2D_colz_graphics()
                    h_2D_res_amp.SetStats(0)
                    h_2D_res_amp.SetYTitle('y [mm]')
                    h_2D_res_amp.SetXTitle('x [mm]')
                    h_2D_res_amp.SetZTitle(var_dT+' [ps]')
                    h_2D_res_amp.Draw('colz')
                    rt.gStyle.SetPaintTextFormat(".0f");
                    h_2D_res_amp.DrawCopy('TEXT SAME')

                    median = np.percentile(ResAmp[:, 0], 50)
                    width = np.percentile(ResAmp[:, 0], 90) - np.percentile(ResAmp[:, 0], 10)
                    sel = np.logical_and(ResAmp[:, 0] < median+2*width, ResAmp[:, 0] > median-2*width)
                    if np.sum(np.logical_not(sel)):
                        print 'Discarting outlayers'
                        print ResAmp[np.logical_not(sel),0]
                    ResAmp = ResAmp[sel]
                    h_2D_res_amp.GetZaxis().SetRangeUser(0.9*np.min(ResAmp[:,0]), 1.1*np.max(ResAmp[:,0]))
                    avg_time_res = np.average(ResAmp[:,0], weights=1./np.square(ResAmp[:,1]))
                    d_avg_time_res = 1./np.sqrt(np.sum(1./np.square(ResAmp[:,1])))
                    results.append([avg_time_res, d_avg_time_res, v_time, v_ref, True])
                    ln = '{:.2f}  {:.2f}  {}  {}  {}\n'.format(avg_time_res, d_avg_time_res, v_time, v_ref, True)
                    file_results.write(ln)
                    if best_result[N_bar].dT[0] < 0 or  best_result[N_bar].dT[0] > avg_time_res:
                        best_result[N_bar].dT = [avg_time_res, d_avg_time_res]
                        best_result[N_bar].var = [v_time, v_ref]
                        best_result[N_bar].AmpCorr = True
                    msg = 'Avg resolution after amp correction: {:.1f} +/- {:.1f} ps'.format(avg_time_res, d_avg_time_res)
                    print msg
                    l = rt.TLatex()
                    l.SetTextSize(0.04);
                    l.DrawLatexNDC(0.15, 0.89, msg.replace('+/-', '#pm'))
                    canvas['c_'+tag].Update()
                    canvas['c_'+tag].SaveAs(out_dir + '/TimeRes2DAmp_bar{:02d}'.format(N_bar)+figform)

        file_results.close()

    print '\n\n======================= Summary =============================='
    table =  PrettyTable(['Bar', 'Best Resolution [ps]', 'Var ref', 'Var timr', 'Amp corrected'])
    for k, res in best_result.iteritems():
        row = [str(k), '{:.2f} +/- {:.2f}'.format(res.dT[0],res.dT[1])]
        row += [res.var[1], res.var[0], 'Yes' if res.AmpCorr else 'No']
        table.add_row(row)

    print table
    table_txt = table.get_string()
    with open(headout_dir+'/SummaryTable.txt','w') as file:
        file.write(table_txt)
