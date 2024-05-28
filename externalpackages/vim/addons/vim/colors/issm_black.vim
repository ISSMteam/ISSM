" issm_black color scheme
" table http://www.calmar.ws/vim/256-xterm-24bit-rgb-color-chart.html

"Set up coloring
hi clear
if exists("syntax_on")
  syntax reset
endif
syntax on
let colors_name = "issm_black"

"preferred colors:
"  0  : black
"  9  : red
" 15  : white
" 69  : blue
" 87  : cyan
" 129 : dark purple
" 202 : orange
" 214 : light orange
" 227 : yellow
" 234 : very dark gray
" 241 : dark gray
" 244 : gray
" 255 : light gray

"                forground    background   style (reverse,bold,..)
" General colors
hi Normal        ctermfg=NONE ctermbg=NONE cterm=NONE
hi NonText       ctermfg=129  ctermbg=NONE cterm=NONE
hi Cursor        ctermfg=NONE ctermbg=NONE cterm=reverse
hi LineNr        ctermfg=15   ctermbg=241  cterm=NONE
hi VertSplit     ctermfg=241  ctermbg=241  cterm=NONE
hi StatusLine    ctermfg=241  ctermbg=87   cterm=NONE
hi StatusLineNC  ctermfg=87   ctermbg=241  cterm=NONE
hi Folded        ctermfg=141  ctermbg=236  cterm=NONE
hi Title         ctermfg=NONE ctermbg=NONE cterm=NONE
hi Visual        ctermfg=NONE ctermbg=NONE cterm=reverse
hi SpecialKey    ctermfg=NONE ctermbg=NONE cterm=NONE
hi WildMenu      ctermfg=0    ctermbg=227  cterm=NONE
hi PmenuSbar     ctermfg=0    ctermbg=129  cterm=NONE
hi Error         ctermfg=15   ctermbg=129  cterm=NONE
hi ErrorMsg      ctermfg=15   ctermbg=129  cterm=NONE
hi WarningMsg    ctermfg=15   ctermbg=129  cterm=NONE

" Message displayed in lower left, such as --INSERT--
hi ModeMsg       ctermfg=0  ctermbg=87   cterm=BOLD
if version >= 700 " Vim 7.x specific colors
  hi CursorLine   ctermfg=NONE ctermbg=NONE cterm=BOLD
  hi CursorColumn ctermfg=NONE ctermbg=NONE cterm=BOLD
  hi MatchParen   ctermfg=87   ctermbg=241  cterm=reverse "matching parenthesis
  hi Pmenu        ctermfg=0    ctermbg=15   cterm=NONE    "auto completion panel
  hi PmenuSel     ctermfg=255  ctermbg=241  cterm=NONE
  hi Search       ctermfg=0    ctermbg=220  cterm=NONE
endif

" Syntax highlighting
hi Comment      ctermfg=244  ctermbg=NONE cterm=NONE
hi String       ctermfg=120  ctermbg=NONE cterm=NONE
hi Number       ctermfg=201  ctermbg=NONE cterm=NONE
hi Keyword      ctermfg=9    ctermbg=NONE cterm=NONE  " matlab function
hi PreProc      ctermfg=9    ctermbg=NONE cterm=NONE  " def undef include
hi Conditional  ctermfg=220  ctermbg=NONE cterm=NONE  " if else end
hi Todo         ctermfg=214  ctermbg=NONE cterm=NONE
hi Constant     ctermfg=196  ctermbg=NONE cterm=NONE
hi Identifier   ctermfg=9    ctermbg=NONE cterm=NONE
hi Function     ctermfg=69   ctermbg=NONE cterm=NONE "functions 69 = pastel blue
hi Type         ctermfg=87   ctermbg=NONE cterm=NONE "cterm matlab global
hi Statement    ctermfg=69   ctermbg=NONE cterm=NONE "cd ls sed mv
hi Special      ctermfg=202  ctermbg=NONE cterm=NONE " matlab '...'
hi Delimiter    ctermfg=NONE ctermbg=NONE cterm=NONE " [ ]
hi Operator     ctermfg=202  ctermbg=NONE cterm=NONE " == &
hi Directory    ctermfg=87   ctermbg=NONE cterm=NONE " == & 

"Specific for diff
hi DiffAdd      cterm=none ctermfg=0 ctermbg=119
hi DiffChange   cterm=none ctermfg=0 ctermbg=228
hi DiffText     cterm=none ctermfg=0 ctermbg=178
hi DiffDelete   cterm=none ctermfg=0 ctermbg=197
hi diffLine     cterm=bold ctermfg=241
hi diffOldLine  cterm=none ctermfg=241
hi diffOldFile  cterm=none ctermfg=241
hi diffNewFile  cterm=none ctermfg=241
hi diffAdded    cterm=none
hi diffRemoved  cterm=none ctermfg=9
hi diffChanged  cterm=none ctermfg=69
