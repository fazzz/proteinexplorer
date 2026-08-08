# PyMOL coloring script (4 group(s))
color gray80, all
select grp_1, 1a8o and ((chain A and resi 209) or (chain A and resi 187) or (chain A and resi 206) or (chain A and resi 208) or (chain A and resi 202) or (chain A and resi 203) or (chain A and resi 214) or (chain A and resi 207) or (chain A and resi 191))
color red, grp_1
# grp_1 = pocket_1
select grp_2, 1a8o and ((chain A and resi 195) or (chain A and resi 197) or (chain A and resi 198) or (chain A and resi 218) or (chain A and resi 155) or (chain A and resi 159) or (chain A and resi 158) or (chain A and resi 161) or (chain A and resi 157) or (chain A and resi 160))
color blue, grp_2
# grp_2 = pocket_2
select grp_3, 1a8o and ((chain A and resi 209) or (chain A and resi 183) or (chain A and resi 187) or (chain A and resi 186) or (chain A and resi 210))
color green, grp_3
# grp_3 = pocket_3
select grp_4, 1a8o and ((chain A and resi 183) or (chain A and resi 179) or (chain A and resi 211) or (chain A and resi 182) or (chain A and resi 186) or (chain A and resi 169))
color yellow, grp_4
# grp_4 = pocket_4
