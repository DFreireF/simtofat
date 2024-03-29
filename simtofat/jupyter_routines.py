from iqtools import *
from lmfit import *
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import os
import matplotlib.pyplot as plt


def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.                        
    return [file[:-1] for file in open(master_filename).readlines()]

def read_and_cut_in_frecuency(filename, lframes, time, skip, xcen, xspan, method = None):
    #get iq object
    iq = get_iq_object(filename)
    #select range to read
    nframes = int(time * iq.fs / lframes)
    sframes = int(skip * iq.fs / lframes)
    #read
    iq.read(nframes = nframes, lframes = lframes, sframes = sframes)
    if method: iq.method = method
    else: iq.method = 'mtm' #'fft', 'mtm', 'welch'
    #create spectrogram
    xx, yy, zz = iq.get_power_spectrogram(nframes, lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
    #cut spectrogram in frecuency
    nxx, nyy, nzz = get_cut_spectrogram(xx, yy, zz, xcen = xcen, xspan = xspan)
    #return array with cutted spectrogram
    return nxx, nyy, nzz

def read_and_get_spectrogram(filename, lframes, time, skip):
    #get iq object
    iq = get_iq_object(filename)
    #select range to read
    nframes = int(time * iq.fs / lframes)
    sframes = int(skip * iq.fs / lframes)
    #read
    iq.read(nframes = nframes, lframes = lframes, sframes = sframes)
    iq.method = 'mtm' #'fft', 'mtm', 'welch'
    #create spectrogram
    xx, yy, zz = iq.get_power_spectrogram(nframes, lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
    return xx, yy, zz

def get_tiq_time(filename):
    #get iq object
    iq = get_iq_object(filename)
    return iq.date_time

def read_and_get_averaged_spectrogram(filename, lframes, time, skip):
    xx, yy, zz = read_and_get_spectrogram(filename, lframes, time, skip)
    axx, _ , azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
    return axx, azz

def get_window_averaged(filename, lframes, time, skip, xcen, xspan, plot = False):
    xx, yy, zz = read_and_cut_in_frecuency(filename, lframes, time, skip, xcen, xspan)
    axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
    if plot:
        plot_spectrum(axx[0,:], azz[0,:], dbm=True)
        plt.show()
    return axx, ayy, azz

def generate_2Dpngs_filelist(filelist, out = '/u/dfreiref/72Br/72Br/2Dfigs/',
                            lframes = 2048, time = 20, skip = 0, cen = 650, span = 5e3):
    file_list = read_masterfile(filelist)

    for index, file in enumerate(file_list):
        xx, yy, zz = read_and_cut_in_frecuency(file, lframes, time, skip, cen, span)
        plot_spectrogram(xx, yy, zz, title = file)
        filename = os.path.basename(file)
        plt.savefig(f'{out}{index}_{filename}.png', dpi = 300, format = 'png')
        plt.close()
        
def generate_1Dpngs_filelist(filelist, out = '/u/dfreiref/72Br/72Br/2Dfigs/',
                            lframes = 2048, time = 20, skip = 0, cen = 650, span = 5e3):
    file_list = read_masterfile(filelist)

    for index, file in enumerate(file_list):
        xx, yy, zz = read_and_cut_in_frecuency(file, lframes, time, skip, cen, span)
        filename = os.path.basename(file)
        axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
        plt.plot(axx[0,:], azz[0,:])
        plt.savefig(f'{out}{index}_{filename}.png', dpi = 300, format = 'png')
        plt.close()
        
def generate_2Dhtmls_filelist(filelist, out = '/u/dfreiref/72Br/72Br/2Dfigs/',
                              lframes = 2048, time = 20, skip = 0, cen = 650, span = 5e3):
    file_list = read_masterfile(filelist)
    for index, file in enumerate(file_list):
        xx, yy, zz = read_and_cut_in_frecuency(file, lframes, time, skip, cen, span)
        fig = plot_interactive_spectrogram(xx, yy, zz, title = file)
        filename = os.path.basename(file)
        fig.write_html(f'{out}{index}_{filename}.html')
        fig.show()
        
def generate_2D1Dshow_filelist(filelist, lframes = 2048, time = 20, skip = 0, cen = 650, span = 5e3):
    file_list = read_masterfile(filelist)

    for index, file in enumerate(file_list):
        xx, yy, zz = read_and_cut_in_frecuency(file, lframes, time, skip, cen, span)
        plot_spectrogram(xx, yy, zz, title = file)
        plt.show()
        axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
        plt.plot(axx[0,:], azz[0,:])
        plt.show()
        
def plot_interactive_spectrogram(xx, yy, zz, title = None, zmin = None, zmax = None, showscale = False, colorscale = 'jet', height = 500, width = 900):
    
    fig = go.Figure(data = go.Heatmap(z = zz, x = xx[0,:], y = yy[:,0],
                                      zmax = zmax, zmin = zmin,
                                      colorscale = colorscale, showscale = showscale))
    
    fig.update_layout(title = title, height = height, width = width,
                      
        xaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text="Frequency (kHz)",
            title_font_size = 19,
            tickfont_size = 15,
            showgrid = True,
            nticks = 10),
                      
        yaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text = "Time (s)",
            title_font_size = 19,
            tickfont_size = 15,
            showgrid = False,
            nticks = 10)
    )
    return fig

def plot_interactive_spectrogram_img(xx, yy, zz, title = 'Spectrogram', dbm = False):
    
    frequency_kHz = [f'{x:0.3f}' for x in xx[0,:]/1000]
    time = [f'{y:0.3f}' for y in yy[:,0]]
    if dbm: power = IQBase.get_dbm(y)
    else: power = zz / zz.max()
    
    panda_df = pd.DataFrame(data = power, 
                            index = time, 
                            columns = frequency_kHz)

    fig = px.imshow(panda_df, color_continuous_scale = 'jet', origin = 'lower')
    fig.update_layout(
        title = title,
        coloraxis_showscale=False,
        xaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text="Frequency (kHz)",
            title_font_size = 18,
            tickfont_size = 9,
            showgrid = True,
            nticks = 10,
        ),
        yaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text = "Time (s)",
            title_font_size = 18,
            tickfont_size = 9,
            showgrid = False,
            nticks = 10
        ),
        coloraxis_colorbar = dict(
            orientation = 'h',
        ),
        height = 500,
        width = 900
        )
    return fig

def plot_interactive_spectrum(x, y, dbm = False):
    frequency_kHz = [f'{x:0.3f}' for x in x/1000]
    if dbm: power = IQBase.get_dbm(y)
    else: power = [p/y.max() for p in y]
    fig = px.line(x = frequency_kHz, y = power, markers = True)
    fig.update_layout(
        xaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text="Frequency (kHz)",
            title_font_size = 20,
            tickfont_size = 16,
            showgrid = True,
            nticks = 12,
        ),
        yaxis = dict(
            showline = True,
            mirror = "ticks",
            linewidth = 2,
            ticks = "inside",
            title_text = "Power Spectral Density",
            title_font_size = 20,
            tickfont_size = 16,
            showgrid = True,
            nticks = 10
        )
    )
    return fig

def study_iso_gs_average_power(file, lframes = 2048, time = 20, skip = 0, xcen = 650, span = 2e3,xceni = -2.7e2, 
                               xceng = 20, xcenb = 1e3, xspan = 3e2, ycen = 10, yspan = 20, every = 3):
    filename = os.path.basename(file)
    xx, yy, zz = read_and_cut_in_frecuency(file, lframes, time, skip, xcen, span)
    plot_spectrogram(xx, yy, zz, title = filename)
    plt.show()

    fig = plot_interactive_spectrogram(xx, yy, zz, title = filename)
    fig.show()
    
    axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
    fig2=plot_interactive_spectrum(axx[0,:], azz[0,:])
    fig2.show()
    
    #iso state
    xxi, yyi, zzi = get_cut_spectrogram(xx, yy, zz, xcen = xceni, 
                                        xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxi, yyi, zzi, title = filename)
    plt.show()
    
    #ground state
    xxg, yyg, zzg = get_cut_spectrogram(xx, yy, zz, xcen = xceng, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxg, yyg, zzg, title = filename)
    plt.show()
    
    #background
    xxb, yyb, zzb = get_cut_spectrogram(xx, yy, zz, xcen = xcenb, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxb, yyb, zzb, title = filename)
    plt.show()
    
    timestep = yy[1,0] - yy[0,0]
    power_frame_average(zzi, zzg, zzb, timestep, every = every)

def study_iso_gs_average_power_light(xx, yy, zz, xceni = -2.7e2, xceng = 20, xcenb = 1e3, xspan = 3e2, ycen = None, yspan = None, every = 3):
    #iso state
    xxi, yyi, zzi = get_cut_spectrogram(xx, yy, zz, xcen = xceni, 
                                        xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxi, yyi, zzi)
    plt.show()
    
    #ground state
    xxg, yyg, zzg = get_cut_spectrogram(xx, yy, zz, xcen = xceng, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxg, yyg, zzg,)
    plt.show()
    
    #background
    xxb, yyb, zzb = get_cut_spectrogram(xx, yy, zz, xcen = xcenb, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxb, yyb, zzb)
    plt.show()
    
    timestep = yy[1,0] - yy[0,0]
    power_frame_average(zzi, zzg, zzb, timestep, every = every)
    
def power_frame_average(ziso, zg, zb, timestep, every = 3):
    apfi = [np.average(ziso[i, :]) for i in range(len(ziso[:,0]))]
    apfg = [np.average(zg[i, :]) for i in range(len(zg[:,0]))]
    apfb = [np.average(zb[i, :]) for i in range(len(zb[:,0]))]

    napfi = list()
    napfg = list()
    napfb = list()
    t = list()
    i = 0
    while (i * every + every <= len(apfi)):
        imin = i * every
        imax = i * every + every
        napfi.append(np.average(apfi[imin:imax]))
        napfg.append(np.average(apfg[imin:imax]))
        napfb.append(np.average(apfb[imin:imax]))
        t.append(imax*timestep)
        i = i+1

    data_dict = dict()
    data_dict['Isomer'] = napfi
    data_dict['Ground State'] = napfg
    data_dict['Background'] = napfb
    df = pd.DataFrame(data=data_dict)

    #total = [x + y for x, y in zip(averaged_power_frame_g, averaged_power_frame_iso)]
    fig = px.line(df, x = t, y = ['Isomer', 'Ground State', 'Background'], markers = True)
    #fig.write_html('/u/litv-exp/personal_directories/dfreiref/72Br/spectrum_isoVSgsVStotal2.html')
    fig.show()
    
def correct_shift(xx, yy, zz, size = 7, cooling = False, heating = False, change_ref = False, show = False):
    # f(t,p) ; t(p,f) ; p(t,f)
    freq_average_per_tframe = dict()
    deltas, nzz = [np.array([]) for i in range(0,2)]
    
    # First we look for our reference / pattern 
    maximum = -1
    for freq_frame, _ in enumerate(zz[0, :]):
        power = np.average(zz[6:, freq_frame])
        if power > maximum:
            maximum = power
            heating_freq = xx[0, freq_frame]
            ref_index = freq_frame
    zero_index = np.argmin(np.abs(xx[0, :]))
    
    for time_frame, _ in enumerate(zz[:, 0]):
        freq_average = np.array([])
        j = 0
        while (j + size <= len(zz[0, :])):
            imin = j
            imax = j + size
            freq_average = np.append(freq_average, np.average(zz[time_frame, imin:imax]))
            j = j + 1
        freq_average_per_tframe[f'{time_frame}'] = freq_average

    for key in reversed(freq_average_per_tframe):
        index_max = np.argmax(freq_average_per_tframe[key][:])
        index_max_f = int((index_max + index_max + size - 1) / 2)
        delta = ref_index - index_max_f
        if change_ref: delta = delta + (zero_index - ref_index)
        nzz = np.append(nzz, np.roll(zz[int(key), :], delta))
        deltas = np.append(deltas, delta)
    nzz = np.reshape(nzz, np.shape(zz))
    
    if cooling: xx = xx + xx[0, int(-deltas[-1] + ref_index)] 
    elif heating: xx = xx + heating_freq 

    deltax = abs(xx[0,0] - xx[0,1])
    deltay = abs(yy[0,0] - yy[1,0])

    shift = {'Time (s)' : np.array([i * deltay for i, _ in enumerate(deltas)])}
    shift['Deltas'] = np.array([delta for delta in deltas])

    a, b, c = fit_decay(shift['Time (s)'], shift['Deltas'])
    fitted_shift = np.array([decay_curve(time, a, b, c) for time in shift['Time (s)']])
    
    distance = [np.abs(fitted_shift[i] - delta) for i, delta in enumerate(shift['Deltas'])]
    for i, dist in enumerate(distance):
        if dist > np.abs(2 * fitted_shift[i]):
            nzz[i] = np.roll(nzz[i], int(dist))
            deltas[i] = int(fitted_shift[i])
    deltas = deltas[::-1]

    shift['Frequency (Hz)'] = np.array([delta * deltax for delta in deltas])
    if show:
        df = pd.DataFrame(shift)
        figi = px.line(df, x = 'Time (s)', y = 'Deltas', markers = True, title = 'Before fit correction')
        figi.show()
        figf = px.line(df, x = 'Time (s)', y = 'Frequency (Hz)', markers = True, title = 'After fit correction')
        figf.show()
        
    return xx, yy[::-1], nzz, deltas

def decay_curve(x, a, b, c):#for the lifetime calculation
    return a + b * np.exp(-x*c) # be careful with /c, it provokes an np.exp()-> infinity overflow

def fit_decay(x, y, seed_a = 0, seed_b = 1, seed_c = 1, fit_report = False):
    gmodel = Model(decay_curve)
    result = gmodel.fit(y, x = x, a = seed_a, b = seed_b, c = seed_c)
    if fit_report:
        fig = px.line(x = x, y = [result.best_fit, y], markers = True)
        fig.show()
        print(result.fit_report())
        
    return result.params['a'].value, result.params['b'].value, result.params['c'].value

def basic_visualization(filename, lframes, time, skip, fcen, fspan):
    xx, yy, zz = read_and_cut_in_frecuency(filename, lframes, time, skip, fcen, fspan)
    axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
    file = os.path.basename(filename)
    plot_spectrogram(xx, yy, zz, title = file)
    plt.show()
    fig = plot_interactive_spectrum(axx[0,:], azz[0,:])
    fig.show()

def power_frame_average_half(ziso, zg, zb, zh1, zh2, timestep, every = 1):
    apfi = [np.average(ziso[i, :]) for i in range(len(ziso[:,0]))]
    apfg = [np.average(zg[i, :]) for i in range(len(zg[:,0]))]
    apfb = [np.average(zb[i, :]) for i in range(len(zb[:,0]))]
    
    apfih1 = [np.average(zh1[i, :]) for i in range(len(zh1[:,0]))]
    apfih2 = [np.average(zh2[i, :]) for i in range(len(zh2[:,0]))]

    napfi = list()
    napfih2 = list()
    napfih1 = list()
    napfg = list()
    napfb = list()
    t = list()
    i = 0
    
    while (i * every + every <= len(apfi)):
        imin = i * every
        imax = i * every + every
        napfi.append(np.average(apfi[imin:imax]))
        
        napfih1.append(np.average(apfih1[imin:imax]))
        napfih2.append(np.average(apfih2[imin:imax]))
        
        napfg.append(np.average(apfg[imin:imax]))
        napfb.append(np.average(apfb[imin:imax]))
        t.append(imax*timestep)
        i = i+1

    data_dict = dict()
    data_dict['Isomer'] = napfi
    data_dict['Isomer_half1'] = napfih1
    data_dict['Isomer_half2'] = napfih2
    data_dict['Ground State'] = napfg
    data_dict['Background'] = napfb
    df = pd.DataFrame(data = data_dict)

    #total = [x + y for x, y in zip(averaged_power_frame_g, averaged_power_frame_iso)]
    fig = px.line(df, x = t, y = ['Isomer', 'Ground State', 'Background', 'Isomer_half1', 'Isomer_half2'], markers = True)
    #fig.write_html('/u/litv-exp/personal_directories/dfreiref/72Br/spectrum_isoVSgsVStotal2.html')
    fig.show()
    return napfi, napfih1, napfih2
    
def study_iso_gs_average_power_half(xx, yy, zz, xceni = -2.7e2, xceng = 20, xcenb = 1e3, xspan = 3e2, ycen = None, yspan = None, every = 1):
    #iso state
    xxi, yyi, zzi = get_cut_spectrogram(xx, yy, zz, xcen = xceni, 
                                        xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxi, yyi, zzi)
    plt.show()
    
    cenh1 = xceni - xspan / 4
    spanh1 = xspan / 2
    cenh2 = xceni + xspan / 4
    spanh2 = xspan / 2
    xxh1, yyh1, zzh1 = get_cut_spectrogram(xx, yy, zz, xcen = cenh1, 
                                        xspan = spanh1, ycen = ycen , yspan = yspan)
    xxh2, yyh2, zzh2 = get_cut_spectrogram(xx, yy, zz, xcen = cenh2, 
                                        xspan = spanh2, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxh1, yyh1, zzh1)
    plt.show()
    plot_spectrogram(xxh2, yyh2, zzh2 )
    plt.show()
    
    #ground state
    xxg, yyg, zzg = get_cut_spectrogram(xx, yy, zz, xcen = xceng, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxg, yyg, zzg,)
    plt.show()
    
    #background
    xxb, yyb, zzb = get_cut_spectrogram(xx, yy, zz, xcen = xcenb, 
                                           xspan = xspan, ycen = ycen , yspan = yspan)
    plot_spectrogram(xxb, yyb, zzb)
    plt.show()
    
    timestep = yy[1,0] - yy[0,0]
    napfi, napfih1, napfih2 = power_frame_average_half(zzi, zzg, zzb, zzh1, zzh2, timestep, every = every)
    
    return napfi, napfih1, napfih2

def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]

def write_spectrum_to_csv(freq, power, filename, center = 0, out = None):
    
    concat_data = np.concatenate((freq, power, IQBase.get_dbm(power)))
    final_data = np.reshape(concat_data, (3, -1)).T
    if out:
        filename = os.path.basename(filename)
    file_name = f'{filename}.csv'
    if out: file_name = os.path.join(out, file_name)
    print(f'created file: {file_name}')
    np.savetxt(file_name, final_data, header =
               f'Delta f [Hz] @ {center} [Hz]|Power [W]|Power [dBm]', delimiter = '|')

def write_tgraph(x,y,file = 'graph'):
    from ROOT import TGraph, TFile
    gra = TGraph(len(np.array(x)),np.array(x),np.array(y))
    ffw = TFile(file+'.root', 'RECREATE')
    gra.Write()
    ffw.Close()
    
def write_etgraph(x,y,ey, ex = None, file = 'egraph'):
    from ROOT import TGraphErrors, TFile
    if ex == None: ex = np.zeros_like(ey)
    egra = TGraphErrors(len(np.array(x)),np.array(x),np.array(y), np.array(ex), np.array(ey))
    ffw = TFile(file+'.root', 'RECREATE')
    egra.Write()
    ffw.Close()
    
def write_spectrum_to_root(xx,yy,zz, filename = 'spectrum'):
    #nxx4w,nyy4w[::-1],total_power4w[::-1]
    from ROOT import TFile
    histo = get_root_th2d(xx,yy,zz)
    file = TFile(filename + '.root', 'RECREATE')
    histo.Write()
    file.Close()
