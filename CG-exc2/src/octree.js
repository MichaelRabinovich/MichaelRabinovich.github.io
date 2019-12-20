var Octree = function (bb, depth) {
    this.bounding = bb; this.data = null;
    this.children = new Array(); this.depth = depth;
	
	this.MAX_OBJ_IN_CELL = 3;
}

Octree.prototype.getIntersectionObjects = function (pos, ray_d) {
    if (this.isLeaf()) {
        if (this.data !== null && this.depth <= 0) return this.data;
        if (this.data !== null) return [this.data];
        return [];
    }
    var objects = new Array();
    for (var i = 0; i < 8; i++) {
        if (this.children[i].bounding.intersects(pos, ray_d)) {
            objects = objects.concat(this.children[i].getIntersectionObjects(pos, ray_d));
        }
    }
    return objects;
}

Octree.prototype.insertObject = function (object) {
    if (this.depth <= 0) {
        if (this.data === null) {
            this.data = new Array();
        }
        this.data.push (object);
        return;
    }

    // The node is a leaf (no children/not split) and has no data assigned
    if (this.isLeaf() && this.data === null) {
        // no data currently assigned and no children,
        // assign this data point to this leaf node
        this.data = object; return;
    }
	
	// The node is a leaf data assigned (but less then max_obj in cell different objects)
    if (this.isLeaf() && this.data.length < this.MAX_OBJ_IN_CELL) {
        // no data currently assigned and no children,
        // assign this data point to this leaf node
        this.data.push(object); return;
    }

    // The node is a leaf but it already has too much objects assigned
    if (this.isLeaf() && this.data !== null) {
        // Split this node into eight children, and then re-insert the old point and
        // our new point into the new children.
		// This could happen several times during insert if these objects are really close
        // to each other.

        this.makeChildren();
        var old_data = this.data;
        this.data = null;
        this.insertObject(old_data);
    }

    // The node is an interior node to the tree (has 8 children)
    if (this.isLeaf() === false) {
        // find out which of the eight children the object collides with then make a recursive call to
        // insert into that child.
        var objBounding = BoundingBox.getObjectBoundingBox(object);

        for (var i = 0; i < this.children.length; i++) {
            if (this.children[i].bounding.collides(objBounding)) {
                this.children[i].insertObject(object);
            }
        }
        return;
    }
}

Octree.prototype.makeChildren = function () {
    for (var x = 0; x <= 1; x++) {
        for (var y = 0; y <= 1; y++) {
            for (var z = 0; z <= 1; z++) {

                var new_b = new BoundingBox( this.bounding.x_max - (1-x)*this.bounding.x_width/2,
											this.bounding.x_min +     x*this.bounding.x_width/2,
											this.bounding.y_max - (1-y)*this.bounding.y_width/2,
											this.bounding.y_min +     y*this.bounding.y_width/2,
											this.bounding.z_max - (1-z)*this.bounding.z_width/2,
											this.bounding.z_min +     z*this.bounding.z_width/2
                )
                this.children.push (new Octree (new_b, this.depth-1));
            }
        }
    }
}

Octree.prototype.isLeaf = function () {
    return this.children.length === 0;
}
