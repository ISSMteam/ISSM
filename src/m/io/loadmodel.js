function loadmodel(modelstring) {

	var md=JSONfn.parse(decodeURI(modelstring));
	md.fix();
	return md;
}
