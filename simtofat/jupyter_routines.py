from iqtools import *
import plotly.express as px
import pandas as pd
import os
import matplotlib.pyplot as plt


def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.                        
    return [file[:-1] for file in open(master_filename).readlines()]

def read_and_cut_in_frecuency(filename, lframes, time, skip, xcen, xspan):
    #get iq object
    iq = get_iq_object(filename)
    #select range to read
    nframes = int(time * iq.fs / lframes)
    sframes = int(skip * iq.fs / lframes)
    #read
    iq.read(nframes = nframes, lframes = lframes, sframes = sframes)
    iq.method = 'mtm' #'fft', 'mtm', 'welch'
    #create spectrogram
    xx, yy, zz = iq.get_spectrogram(nframes, lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
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
    xx, yy, zz = iq.get_spectrogram(nframes, lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
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
        
def plot_interactive_spectrogram(xx, yy, zz, title = 'Spectrogram'):
    frequency_kHz = [f'{x:0.3f}' for x in xx[0,:]/1000]
    time = [f'{y:0.3f}' for y in yy[:,0]]
    norm_power = zz / zz.max()
    
    panda_df = pd.DataFrame(data = norm_power, 
                            index = time, 
                            columns = frequency_kHz)

    fig = px.imshow(panda_df, color_continuous_scale = 'jet', origin = 'lower')
    fig.update_layout(
        title = title,
        #font_family="Avenir",
        #hoverlabel_font_family="Avenir",
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
        height = 1080,
        width = 1980
        )
    return fig

def plot_interactive_spectrum(x, y):
    frequency_kHz = [f'{x:0.3f}' for x in x/1000]
    norm_power = [p/y.max() for p in y]
    fig = px.line(x = frequency_kHz, y = norm_power, markers = True)
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
    
def correct_shift(zz, every = 7, heating = False):
    # xx, yy, zz -> Spectrogram variables
    # f(t,p) ; t(p,f) ; p(t,f)
    freq_average_per_tframe = dict()
    ref = -1
    deltas, nzz = [np.array([]) for i in range(0,2)]
    
    for time_frame, _ in enumerate((zz[:,0])):
        freq_average = np.array([])
        i = 0
        while (i + every <= len(zz[0,:])):
            imin = i
            imax = i + every
            freq_average = np.append(freq_average, np.average(zz[time_frame, imin:imax]))
            freq_average_per_tframe[f'{time_frame}'] = freq_average
            i = i + 1
            
    if heating: freq_average_per_frame = reversed(freq_average_per_frame)
    for key in freq_average_per_tframe:
        max = freq_average_per_tframe[key][:].max()
        min = freq_average_per_tframe[key][:].min()
        if max >= 1.4 * min:
            index_max = np.argmax(freq_average_per_tframe[key][:])
            index_max_f = int((index_max + index_max + every - 1) / 2) # This has to be integer, if every is even (like 5,7,9)
            if ref == -1: ref = index_max_f
            else: delta = ref - index_max_f
            nzz = np.append(nzz, np.roll(zz[int(key), :], delta))
            deltas = np.append(deltas, delta)
        else: nzz = np.append(nzz, zz[int(key), :])
    nzz = np.reshape(nzz, np.shape(zz))
    if heating: return reversed(nzz), reversed(deltas)
    else: return nzz, deltas

def basic_visualization(filename, lframes, time, skip, fcen, fspan):
    xx, yy, zz = read_and_cut_in_frecuency(filename, lframes, time, skip, fcen, fspan)
    axx, ayy, azz = get_averaged_spectrogram(xx, yy, zz, len(xx[:,0]))
    file = os.path.basename(filename)
    plot_spectrogram(xx, yy, zz, title = file)
    plt.show()
    fig = plot_interactive_spectrum(axx[0,:], azz[0,:])
    fig.show()
