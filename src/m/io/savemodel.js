function savemodel(md) {

	var string=encodeURI(JSONfn.stringify(md));

	var url='data:text/json:charset=utf8,' + encodeURIComponent(string);
	window.open(url, '_blank');
	window.focus();
}
