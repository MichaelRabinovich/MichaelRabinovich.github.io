var BoundingBox = function (_x_max, _x_min, _y_max, _y_min, _z_max, _z_min) {
    this.x_max = _x_max;
    this.x_min = _x_min;
    this.y_max = _y_max;
    this.y_min = _y_min;
    this.z_max = _z_max;
    this.z_min = _z_min;

    this.x_width = this.x_max - this.x_min;
    this.y_width = this.y_max - this.y_min;
    this.z_width = this.z_max - this.z_min;
}

BoundingBox.prototype.collides = function (bounding) {
    if (bounding.x_min > this.x_max || bounding.x_max < this.x_min) return false;
    if (bounding.y_min > this.y_max || bounding.y_max < this.y_min) return false;
    if (bounding.z_min > this.z_max || bounding.z_max < this.z_min) return false;

    return true;
}

BoundingBox.prototype.intersects = function (pos, ray_d) {
    var t_x1 = (this.x_min - pos.e(1)) / ray_d.e(1);
    var t_y1 = (this.y_min - pos.e(2)) / ray_d.e(2);
    var t_z1 = (this.z_min - pos.e(3)) / ray_d.e(3);
    var t_x2 = (this.x_max - pos.e(1)) / ray_d.e(1);
    var t_y2 = (this.y_max - pos.e(2)) / ray_d.e(2);
    var t_z2 = (this.z_max - pos.e(3)) / ray_d.e(3);

    var t_x_min = Math.min (t_x1, t_x2);
    var t_x_max = Math.max (t_x1, t_x2);
    var t_y_min = Math.min (t_y1, t_y2);
    var t_y_max = Math.max (t_y1, t_y2);
    var t_z_min = Math.min (t_z1, t_z2);
    var t_z_max = Math.max (t_z1, t_z2);

    var t_min = Math.max (t_x_min, t_y_min, t_z_min);
    var t_max = Math.min (t_x_max, t_y_max, t_z_max);

    // if t_min < t_max it intersects
    return (t_min < t_max);
}

BoundingBox.getMultipleObjectsBoundingBox = function (objects) {
    var x_min = Infinity;
    var y_min = Infinity;
    var z_min = Infinity;
    var x_max = -Infinity;
    var y_max = -Infinity;
    var z_max = -Infinity;

    for (var i = 0; i < objects.length; i++) {
        var bounding = BoundingBox.getObjectBoundingBox(objects[i]);
		x_min = Math.min(x_min, bounding.x_min);
		y_min = Math.min(y_min, bounding.y_min);
		z_min = Math.min(z_min, bounding.z_min);
        
		x_max = Math.max(x_max,bounding.x_max);
		y_max = Math.max(y_max,bounding.y_max);
		z_max = Math.max(z_max,bounding.z_max);
        
    }
    return new BoundingBox(x_max,x_min,y_max,y_min,z_max,z_min);
}

BoundingBox.getObjectBoundingBox = function (obj) {
	if (obj.Type == "Triangle" ) {
		var x_min = Math.min(obj.v1.e(1),obj.v2.e(1),obj.v3.e(1)); var x_max = Math.max(obj.v1.e(1),obj.v2.e(1),obj.v3.e(1));
		var y_min = Math.min(obj.v1.e(2),obj.v2.e(2),obj.v3.e(2)); var y_max = Math.max(obj.v1.e(2),obj.v2.e(2),obj.v3.e(2));
		var z_min = Math.min(obj.v1.e(3),obj.v2.e(3),obj.v3.e(3)); var z_max = Math.max(obj.v1.e(3),obj.v2.e(3),obj.v3.e(3));
		return new BoundingBox(x_max,x_min,y_max,y_min,z_max,z_min);
	} else if (obj.Type == "Sphere") {
		return new BoundingBox(obj.Center.e(1) + obj.Radius, obj.Center.e(1) - obj.Radius,
							   obj.Center.e(2) + obj.Radius, obj.Center.e(2) - obj.Radius,
							   obj.Center.e(3) + obj.Radius, obj.Center.e(3) - obj.Radius)
	}
	// Should never get here
}