var Triangle = function (v1, v2, v3, ambientColor, diffuseColor,specularColor, specularExponent) {
	this.Type = "Triangle";
	
    this.v1 = v1; this.n1 = null;
    this.v2 = v2; this.n2 = null;
    this.v3 = v3; this.n3 = null;
	
	//Find vectors for two edges sharing V1
	this.e1 = this.v2.subtract(v1); this.e2 = this.v3.subtract(v1);
	this.normal = this.e1.cross(this.e2).toUnitVector();
	
	this.AmbientColor = ambientColor;
	this.DiffuseColor = diffuseColor;
	this.SpecularColor = specularColor;
	this.SpecularExponent = specularExponent;
};


Triangle.prototype.getRayPtNormal = function(intersectionPoint) {
    var a1 = new Triangle(intersectionPoint, this.v2, this.v3).get_area();
    var a2 = new Triangle(intersectionPoint, this.v1, this.v3).get_area();
    var a3 = new Triangle(intersectionPoint, this.v1, this.v2).get_area();

    var normal = $V(0,0,0);
    normal = normal.add (this.n1.multiply(a1));
    normal = normal.add (this.n2.multiply(a2));
    normal = normal.add (this.n3.multiply(a3));

    return normal.toUnitVector();
}

Triangle.prototype.get_area = function () {
	var AB = this.v2.subtract(this.v1);
	var AC = this.v3.subtract(this.v1);
	var cr = AB.cross(AC);
	return Math.sqrt (cr.e(1)*cr.e(1) + cr.e(2)*cr.e(2) + cr.e(3)*cr.e(3)) / 2;
}

Triangle.prototype.intersect = function (pos, ray) {

	EPSILON = 0.00001;
	v1_v2 = this.v2.subtract(this.v1);
	v1_v3 = this.v1.subtract(this.v3);
	n = this.normal;
	nDotRay = n.dot(ray);
	if (Math.abs(nDotRay) < EPSILON) return null;
	var d = n.dot(this.v1);
	var t = (-n.dot(pos) + d )/nDotRay;
	Phit = pos.add(ray.multiply(t));
    // inside-out test
  
	// inside-out test edge0
	var v1p = Phit.subtract(this.v1);
	var v = n.dot(v1_v2.cross(v1p));
	if (v < 0) return null; // P outside triangle
	// inside-out test edge1
	var v2p = Phit.subtract(this.v2);
	var v2v3 = this.v3.subtract(this.v2);
    var w = n.dot(v2v3.cross(v2p));
	if (w < 0) return null; // P outside triangle
 
	// inside-out test edge2
	var v3p = Phit.subtract(this.v3);
	var v3v1 = this.v1.subtract(this.v3);
	var u = n.dot(v3v1.cross(v3p));
	if (u < 0) return null; // P outside triangle
	
	return t;
}