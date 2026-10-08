# HIJING 0-10% background images for the overlay samples

Source: yeonjugo's `/sphenix/user/yeonjugo/jetML/create_calo_images/macro/type4_run19_hijingAll_noNoise_0to10/histograms/cent0/`
(`type4_run19_hijingAll_noNoise_images_cent0_file<N>.root`, N = 0..49999): HIJING Au+Au 0-10%,
no calorimeter noise, one 24x64 eta-phi TH2D per event, made with the same `draw_event_v2.C`
image algorithm as our images.

`type4_run19_hijingAll_noNoise_0to10_cent0.counts`: `<file N> <number of images>` for every file
(counted 2026-10-05): 985,791 images, 3 empty files, none unreadable.

HIJING event e (e = 0, 1, ...) is image number p (in key order) of file N, where files are taken in
order N = 0, 1, ... and e = (images in files before N) + p. `macros/overlay_hijing.C` adds HIJING
event e to JEWEL event e of a sample (JEWEL events numbered in pass-1 segment order). Without
reuse, the overlay samples hold the first 985,791 JEWEL events.
