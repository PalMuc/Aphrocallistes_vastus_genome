hide everything
show cartoon
set_color colordefault, [0.75,0.75,0.58]
color colordefault, all
set_color blue0, [0.63,0.63,0.63]
set_color blue50, [0.5,0.58,0.68]
set_color blue60, [0.42,0.55,0.71]
set_color blue70, [0.35,0.52,0.73]
set_color blue80, [0.28,0.49,0.76]
set_color blue90, [0.2,0.46,0.8]
set_color blue95, [0.12,0.43,0.83]
set_color blue98, [0.0,0.38,0.87]
set_color blue100, [0.6,0,1]
color blue0, chain A
select 50pct_grp_2_A, (chain A & resi 37,65,137,138,143,147,170,207,208,233,272)
color blue50, 50pct_grp_2_A
select 60pct_grp_3_A, (chain A & resi 77,106,122,169,172,187,224,268,269)
color blue60, 60pct_grp_3_A
select 70pct_grp_4_A, (chain A & resi 88,101,115,117,120,148,163,213,222,227,257)
color blue70, 70pct_grp_4_A
select 80pct_grp_5_A, (chain A & resi 55,61,119,121,132,146,167,220,225,267)
color blue80, 80pct_grp_5_A
select 90pct_grp_6_A, (chain A & resi 39,63,125,134,216,223,229,265)
color blue90, 90pct_grp_6_A
select 95pct_grp_7_A, (chain A & resi 45,124,221)
color blue95, 95pct_grp_7_A
select 98pct_grp_8_A, (chain A & resi 123,219,260,262,270)
color blue98, 98pct_grp_8_A
select 100pct_grp_9_A, (chain A & resi 58,59,60,131,133,144,218,231)
color blue100, 100pct_grp_9_A
show sticks, 100pct_grp_9_A
select sele, (chain A & resi 1-30)
hide everything, sele
select active_site_H, (chain A & resi 121,123,146)
color paleyellow, active_site_H
show sticks, active_site_H
select active_site_T, (chain A & resi 221,222)
color red, active_site_T
show sticks, active_site_T
select active_site_Y, (chain A & resi 95)
color tv_orange, active_site_Y
show sticks, active_site_Y
deselect
bg white
