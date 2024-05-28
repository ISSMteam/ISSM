function! IssmFoldText()

	" Get line {{{
	let line = getline(v:foldstart)
	"let line = line . '--'
	"}}}
	" remove the marker that caused this fold from the display {{{
	let foldmarkers = split(&foldmarker,',')
	let line = substitute(line, '\V' .  foldmarkers[0] . '\%(\d\+\)\?', ' ', '')
	" }}}
	" remove comments that vim knows about {{{
	let comment = split(&commentstring, '%s')
	if comment[0] != ''
		let comment_begin = comment[0]
		let comment_end = ''
		if len(comment) > 1
			let comment_end = comment[1]
		end
		let pattern = '\V' .  comment_begin .  '\s\*' .  comment_end .  '\s\*\$'
		if line =~ pattern
			let line = substitute(line, pattern, ' ', '')
		else
			let line = substitute(line, '.*\V' .  comment_begin, ' ', '')
			if comment_end != ''
				let line = substitute(line, '\V' .  comment_end, ' ', '')
			endif
		endif
	endif
	" }}}
	" remove any remaining leading or trailing whitespace {{{
	"let line = substitute(line, '^\s*\(.\{-}\)\s*$', '\1', '')
	let line = substitute(line, '^\s*%\(.\{-}\)\s*$', '\1', '') "Also remove % in matlab comments
	" }}}
	" align everything, and pad the end of the display with - {{{
	let alignment = &columns - 18 - v:foldlevel
	let line = strpart(printf('%-'.alignment.'s',line),0,alignment)
	"let line = substitute(line, '\%( \)\@<= \%( *$\)\@=', '-', 'g') " ->dashes
	let line = substitute(line, '\%( \)\@<= \%( *$\)\@=',' ', 'g')  " ->white spaces
	" }}}
	" format the line count {{{
	let cnt = printf('%13s','('.(v:foldend - v:foldstart + 1) .' lines) ')
	" }}}
	return '+-'.v:folddashes.' '.line.cnt

endfunction
