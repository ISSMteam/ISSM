function supportedcostfunctions() {
    function range(start, count) {
      return Array.apply(0, Array(count))
        .map(function (element, index) { 
          return index + start;  
      });
    }
    var list = range(101,5).concat(range(201,1).concat(range(501,7)).concat(range(510,1)).concat(range(601,4)));
    return list;
}
