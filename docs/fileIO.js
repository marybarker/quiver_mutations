function saveFile(fmt='globalQP') {
  var data = JSON.stringify({
    "nodes": QPglobalNodes,
    "edges": QPglobalEdges,
    "frozenNodes": QPglobalFrozenNodes,
    "potential": QPglobalPotential
  });

  if (fmt == 'globalTriangulation') {
    var opt = makeTriangulation(TRIglobalEdges, TRIglobalCoords, TRIglobalBoundaryEdges);
    data = JSON.stringify({
      "nodes": TRIglobalCoords,
      "edges": TRIglobalEdges,
      "triangles": opt[0]
    });
  } else if (fmt == 'QPFromTriangulation') {
    data = JSON.stringify(QPFromTriangulation(TRIglobalTriangulation));
  } else {
    data = allUniqueTriangulations(TRIglobalTriangulation, TRIglobalBoundaryEdges);
    const cs = data.coordinates.map(y => y.map(z=>parseFloat(z)));
    data = data.triangulations.map(x => JSON.parse(x));
    data = data.map(t => t.map(e => JSON.parse(e)));
    data = data.map(x => QPFromTriangulation([x, cs]));

    data = JSON.stringify(data)
  }

  const textToBLOB = new Blob([data], {type: "text/plain"});
  const sFileName = "outputs.JSON";

  let newLink = document.createElement("a");
  newLink.download = sFileName;

  if (window.webkitURL != null) {
    newLink.href = window.webkitURL.createObjectURL(textToBLOB);
  }
  else {
    newLink.href = window.URL.createObjectURL(textToBLOB);
    newLink.style.display = "none";
    document.body.appendChild(newLink);
  }
  newLink.click();
}

let uploadFile = () => {
  var files = document.getElementById("jsonupload").files;
  if (files.length <= 0) {
    return false;
  }
  var fr = new FileReader();
    fr.onload = function(e) { 
    var result = JSON.parse(e.target.result);
    updateQPFromJSON(result);
  }
  fr.readAsText(files.item(0));
};
