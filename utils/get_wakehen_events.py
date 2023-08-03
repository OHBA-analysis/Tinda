import mne
from scipy.io import savemat
from glob import glob
import numpy as np
from osl_dynamics.inference.modes import convert_to_mne_raw
import pickle


fnames = np.sort(glob("/well/woolrich/projects/wakeman_henson/spring23/src/*/sflip_parc-raw.fif"))
alphas = pickle.load(open("/well/woolrich/users/wbs888/wakeman-henson/alp.pkl", "rb"))

# create an events array with the timings relative to the timeline excluding bad segments
events=[]
for alpha, fname in zip(alphas, fnames):
    raw = mne.io.read_raw(fname)
    raw = convert_to_mne_raw(alpha, raw, n_embeddings=15)
    event = mne.find_events(raw, min_duration=0.005)
    event[:,0] -= raw.first_samp
    
    # create event time series with just the number of the event
    event_ts = np.zeros(raw.get_data().shape[1], dtype=np.int64)
    for i, e in enumerate(event[:,0]):
        event_ts[e] = i+1
        
    # Get time indices excluding bad segments from raw
    _, times = raw.get_data(
        reject_by_annotation="omit", return_times=True)
    indices = raw.time_as_index(times)
    
    # get the event time course without bad segments
    event_ts_crop = [event_ts[j] for j in indices]
    
    # find the events in the event time course that still exist
    new_event_idx = np.where([e!=0 for e in event_ts_crop])[0]
    
    # remove the ones included in bad segments
    keep = [event_ts_crop[e]-1 for e in new_event_idx]
    event = event[keep,:]
    
    # substitute the badsegment-corrected trial onsets
    event[:,0] = new_event_idx
    events.append(event)
    
savemat("/well/woolrich/users/wbs888/wakehan_events_badseg_adapted.mat", {"events":events})
