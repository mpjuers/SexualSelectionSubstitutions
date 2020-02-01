let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 Scripts/Batch/execSnakemake.batch
badd +44 term://.//33588:/usr/local/bin/fish
badd +0 Snakefile
badd +0 clusterConfig.json
badd +0 snakemakeConfig.json
badd +5 ~/Code/Modules/seqtk
badd +1 ~/Code/Modules/module_template.tcl
badd +9 term://.//34516:/usr/local/bin/fish
badd +5 ~/Code/Modules/seqtk/1.3
badd +12 term://.//35987:/usr/local/bin/fish
argglobal
%argdel
$argadd Scripts/Batch/execSnakemake.batch
edit snakemakeConfig.json
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd _ | wincmd |
split
wincmd _ | wincmd |
split
2wincmd k
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
wincmd w
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe '1resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 1resize ' . ((&columns * 139 + 180) / 360)
exe '2resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 2resize ' . ((&columns * 139 + 180) / 360)
exe '3resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 3resize ' . ((&columns * 279 + 180) / 360)
exe '4resize ' . ((&lines * 30 + 48) / 97)
exe 'vert 4resize ' . ((&columns * 279 + 180) / 360)
exe 'vert 5resize ' . ((&columns * 80 + 180) / 360)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
1
normal! zo
3
normal! zo
let s:l = 2 - ((1 * winheight(0) + 15) / 31)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
lcd /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions
wincmd w
argglobal
if bufexists("/Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/clusterConfig.json") | buffer /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/clusterConfig.json | else | edit /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/clusterConfig.json | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
1
normal! zo
let s:l = 12 - ((11 * winheight(0) + 15) / 31)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
12
normal! 0
lcd /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions
wincmd w
argglobal
if bufexists("/Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Snakefile") | buffer /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Snakefile | else | edit /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Snakefile | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 13 - ((8 * winheight(0) + 15) / 31)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
13
normal! 0
lcd /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions
wincmd w
argglobal
if bufexists("/Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Scripts/Batch/execSnakemake.batch") | buffer /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Scripts/Batch/execSnakemake.batch | else | edit /Volumes/Carbonate/gpfs/home/m/p/mpjuers/Carbonate/Projects/SexualSelectionSubstitutions/Scripts/Batch/execSnakemake.batch | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 15) / 30)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
if bufexists("term://.//33588:/usr/local/bin/fish") | buffer term://.//33588:/usr/local/bin/fish | else | edit term://.//33588:/usr/local/bin/fish | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 46 - ((45 * winheight(0) + 47) / 94)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
46
normal! 042|
wincmd w
5wincmd w
exe '1resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 1resize ' . ((&columns * 139 + 180) / 360)
exe '2resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 2resize ' . ((&columns * 139 + 180) / 360)
exe '3resize ' . ((&lines * 31 + 48) / 97)
exe 'vert 3resize ' . ((&columns * 279 + 180) / 360)
exe '4resize ' . ((&lines * 30 + 48) / 97)
exe 'vert 4resize ' . ((&columns * 279 + 180) / 360)
exe 'vert 5resize ' . ((&columns * 80 + 180) / 360)
tabnext 1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToOF
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
