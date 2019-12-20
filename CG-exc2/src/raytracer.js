// global variables

// objects
var LargeRedSphere = {
	Type: "Sphere",
	Center: $V(0,0,0),
	Radius: 2,
	AmbientColor: $V(0.75,0,0),
	DiffuseColor: $V(1,0,0),
	SpecularColor: $V(1,1,1),
	SpecularExponent: 32.0,
	b1_refraction: false,
	texture:null, normal_map:null
};

var SmallBlueSphere = {
	Type: "Sphere",
	Center: $V(1.25,1.25,3),
	Radius: 0.5,
	AmbientColor: $V(0,0,0.75),
	DiffuseColor: $V(0,0,1),
	SpecularColor: $V(0.5,0.5,1),
	SpecularExponent: 16.0,
	texture:null, normal_map:null
};

var SmallBlueSphere2 = {
	Type: "Sphere",
	Center: $V(0.25,1.25,3),
	Radius: 1,
	AmbientColor: $V(0,0,0.75),
	DiffuseColor: $V(0,0,1),
	SpecularColor: $V(0.5,0.5,1),
	SpecularExponent: 16.0,
}

var SphereIntersection = {
	Sphere1:SmallBlueSphere,
	Sphere2:SmallBlueSphere2
};


var HemiSphere_x_z = {
	Type: "Hemisphere",
	Center: $V(0,0,0),
	Radius: 2,
	Plane_Normal: $V(-1,0,1).toUnitVector(),
	Plane_Pt: $V(0,0,0), // just because I'm lazy
	Ext_AmbientColor: $V(0.75,0,0),
	Ext_DiffuseColor: $V(1,0,0),
	Ext_SpecularColor: $V(1,1,1),
	Ext_SpecularExponent: 32.0,
	Int_AmbientColor: $V(0.75,0.75,0),
	Int_DiffuseColor: $V(1,1,0),
	Int_SpecularColor: $V(1,1,1),
	Int_SpecularExponent: 32.0,
}

// probably better to keep it inside the HemiSphere
var YellowPlane = {
	Plane_Normal: $V(-1,0,1).toUnitVector(),
	Plane_Pt: $V(0,0,0), // just because I'm lazy
	AmbientColor: $V(0.75,0.75,0),
	DiffuseColor: $V(1,1,0),
	SpecularColor: $V(1,1,1),
	SpecularExponent: 32.0,
}

// Globals
var G = {
	// Scene
	pic_pixels_w:800, pic_pixels_h:600,
	view_pos:$V(0,0,10),
	view_pose_stereo_l: $V(-0.5,0,10),
	view_pose_stereo_r: $V(0.5,0,10),
	view_direction:$V(0,0,-1),
	up_direction:$V(0,1,0),
	proj_plane_d:1, proj_plane_z:9,
	field_of_view:40, // in degrees
	global_ambient:0.2,
	light_p:$V(10,10,10),
	light_c:$V(1,1,1),
	proj_plane_min_x:0, proj_plane_min_y:0, // uninitialized
	plane_pixel_size:0, // uninitialized
	SceneObjects:[],
	Use_Octree:1, octree:null,
	AREA_LIGHT_RADIUS:1, AREA_LIGHT_SAMPLES_NUM:75,
	EPSILON:0.0001,
};

// 0. set up the scene described in the exercise sheet (this is called before the rendering loop)
	
function loadScene() {
	var proj_plane_h = G.proj_plane_d * get_tan_deg(G.field_of_view/2);
	var proj_plane_w = proj_plane_h * (G.pic_pixels_w/G.pic_pixels_h);
	
	G.proj_plane_min_x = G.view_pos.e(1) -proj_plane_w; G.proj_plane_min_y = G.view_pos.e(2) - proj_plane_h;
	G.plane_pixel_size = 2*proj_plane_w / G.pic_pixels_w;
	
	// move first coords to the center
	G.proj_plane_min_x += G.plane_pixel_size; G.proj_plane_min_y += G.plane_pixel_size;
	
	
	// Create scene objects
	if (ModuleId.C3) {
		 var triMesh = readOBJ('./data/sphere.obj');
		 var r_t_scale = 2;
		 var b_t_scale = 0.5;var b_t_translate = $V(1.25,1.25,3);
		 for (var i = 0; i < triMesh.F.length; i++) {
			var f = triMesh.F[i];
			var r_t = new Triangle(triMesh.V[f.e(1)].multiply(r_t_scale),triMesh.V[f.e(2)].multiply(r_t_scale),triMesh.V[f.e(3)].multiply(r_t_scale), $V(0.75,0,0), $V(1,0,0), $V(1,1,1),32.0);
			r_t.n1 = triMesh.N[f.e(1)];
			r_t.n2 = triMesh.N[f.e(2)];
			r_t.n3 = triMesh.N[f.e(3)];
			
			var b_t = new Triangle(triMesh.V[f.e(1)].multiply(b_t_scale).add(b_t_translate),triMesh.V[f.e(2)].multiply(b_t_scale).add(b_t_translate),triMesh.V[f.e(3)].multiply(b_t_scale).add(b_t_translate), $V(0,0,0.75), $V(0,0,1), $V(0.5,0.5,1),16.0);
			b_t.n1 = triMesh.N[f.e(1)];
			b_t.n2 = triMesh.N[f.e(2)];
			b_t.n3 = triMesh.N[f.e(3)];
			
			G.SceneObjects.push(r_t);
			G.SceneObjects.push(b_t);
		}
	} else if (ModuleId.C1) {
		G.SceneObjects.push(HemiSphere_x_z);
		G.SceneObjects.push(SphereIntersection);
		G.Use_Octree = 0;
	} else if (ModuleId.D1)  {		
		for (var i = 0; i <= 9; i++) {
			for (var j = 0; j <= 9; j++) {
				for (var k = 0; k <= 9; k++) {
					var D1BlueSphere = {
						Type: "Sphere",
						Center: $V(i-4.5,j-4.5,-Math.pow(k,3)),
						Radius: 0.25,
						AmbientColor: $V(0,0,0.75),
						DiffuseColor: $V(0,0,1),
						SpecularColor: $V(0.5,0.5,1),
						SpecularExponent: 16.0,
					};
					G.SceneObjects.push(D1BlueSphere);
				}
			}
		}
	} else {
		G.SceneObjects.push(LargeRedSphere);
		G.SceneObjects.push(SmallBlueSphere);
	}
	
	if (ModuleId.C2) {
		LargeRedSphere.texture = readTGA('./data/Earth.tga'); LargeRedSphere.normal_map = readTGA('./data/EarthNormal.tga');
		SmallBlueSphere.texture = readTGA('./data/Moon.tga'); SmallBlueSphere.normal_map = readTGA('./data/MoonNormal.tga');
	}
	
	if (G.Use_Octree) {
        var depth = Math.floor(Math.log (G.SceneObjects.length)/ Math.log(4));
        G.octree = new Octree(BoundingBox.getMultipleObjectsBoundingBox(G.SceneObjects), depth);
        for (var i = 0; i < G.SceneObjects.length; i++) {
            G.octree.insertObject(G.SceneObjects[i]);
        }
    }
}

function trace(color, pixelX, pixelY) {
	var f_color = trace_normal(color, pixelX, pixelY, G.view_pos);
	
	// 5. set the pixel color into the image buffer using the computed shading (for now set dummy color into the image buffer)
	color.setElements(f_color.e(1),f_color.e(2),f_color.e(3));	
}

function supersample(color, pixelX, pixelY, view_pos) {
	G.pic_pixels_w = G.s_pic_pixels_w; G.pic_pixels_h = G.s_pic_pixels_h; G.plane_pixel_size = G.s_plane_pixel_size;
	var f_color = $V(0,0,0); var x_s_step = 0; var y_s_step = 0;
	for (x_s_step = 0; x_s_step < 4; x_s_step++) {
		for (y_s_step = 0; y_s_step < 4; y_s_step++) {
			f_color.addN(trace_normal(color, (pixelX)*4 + x_s_step, (pixelY)*4 + y_s_step, view_pos));
		}
	}
	return f_color.multiply(1/16.);
}
function trace_normal(color, pixelX, pixelY, view_pos) {
	var f_color = $V(0,0,0);
	// 1. shoot a ray determined from the camera parameters and the pixel position in the image
	pixelY = G.pic_pixels_h - pixelY; // opposite y
	var proj_pixel_pos = $V(G.proj_plane_min_x + (pixelX) * G.plane_pixel_size, G.proj_plane_min_y + (pixelY) * G.plane_pixel_size, G.proj_plane_z);
	var ray_d = proj_pixel_pos.subtract(view_pos); ray_d.toUnitVectorN();

	var ret = shoot_ray_and_calc_phong(view_pos, ray_d); f_color  = ret[0];
	return f_color;
}

function shoot_ray_and_calc_phong(pos, ray_d) {
	var f_color = $V(0,0,0); var d_color_v = $V(0,0,0); var spec_color = $V(0,0,0); var amb_color = $V(0,0,0); var n = null;
	
	// 2. intersect the ray to scene elements and determine the closest one
	var ret = find_intersections(pos, ray_d, G.SceneObjects); intersec_obj = ret[0]; intersec_loc = ret[1]; n = ret[2];
	
	if (intersec_obj != null) {
		ray_d.multiplyN(-1);
		if (ModuleId.C2 == true) {
			uv = get_sphere_uv(n); u = uv[0]; v = uv[1];
			text_color = get_tga_bilinear_rgb_by_x_y(intersec_obj.texture,u,v);
			n = get_sphere_normal_from_normal_map(intersec_obj.normal_map, n, u, v);
			 
			mat_ambient = text_color; mat_diffuse = text_color; mat_spec = text_color;
		} else {
			mat_ambient = intersec_obj.AmbientColor; mat_diffuse = intersec_obj.DiffuseColor; mat_spec = intersec_obj.SpecularColor;
		}
		amb_color = mat_ambient.multiply(G.global_ambient);
		
		// 3. check if the intersection point is illuminated by each light source
		var shadow_ray = G.light_p.subtract(intersec_loc); shadow_ray.toUnitVectorN(); 
		
		if (ModuleId.D2) {
			var light_shadow_scalar = getSoftShadow(ray_d, intersec_loc, G.AREA_LIGHT_RADIUS, G.AREA_LIGHT_SAMPLES_NUM);
		} else {
			var light_shadow_scalar = 1;
			G.light_p.subtract(intersec_loc); shadow_ray.toUnitVectorN(); 
			ret = find_intersections(intersec_loc, shadow_ray, G.SceneObjects); intersec_shadow_obj = ret[0]; intersec_shadow_loc = ret[1];
		}
		
		if (ModuleId.D2 || (ModuleId.D2 == false && intersec_shadow_obj == null) ) {
			// 4. shade the intersection point using the material attributes and the lightings
			
			d_color_v = $V(0,0,0);
			var E = shadow_ray.dot(n);
			d_color_v = mat_diffuse.multiply(E).multiply(light_shadow_scalar);
			
			var reflected_spec_ray = n.multiply(2*(shadow_ray.dot(n))).subtract(shadow_ray); reflected_spec_ray.multiply(light_shadow_scalar).toUnitVectorN();
			spec_color = mat_spec.multiply(Math.pow(reflected_spec_ray.dot(ray_d),intersec_obj.SpecularExponent));
		}
		f_color = d_color_v.add(amb_color).add(spec_color);
	}
	return [f_color, intersec_obj, intersec_loc, n];
}

function get_sphere_uv(n) {
	// (0,0,1)->(1,0), (-1,-1,0)->(0,1), therefore we should find a rotation such as (0,1,1)->(0,0,1) and (-1,1,-1)->(-1,1,-1)
	n2 = $V2([n.e(1),n.e(2), n.e(3)]).toUnitVector(); var orig_pole = $V2([0,1,0]);var orig_mer = $V2([-1,0,0]);
	angle_pole = -Math.acos( orig_pole.dot($V2([0,1,1])) );
	n = n2.rotate(angle_pole, $L($V2([0,0,0]), orig_pole) );
	orig_mer = orig_mer.rotate(angle_pole, $L($V2([0,0,0]), orig_pole) );
	angle_mer = -Math.acos( orig_mer.dot($V2([-1,1,-1])) );
	n = n.rotate(angle_mer, $L($V2([0,0,0]), orig_mer) );
	
	x = 0.5 + Math.atan2(n.e(1),n.e(3)) / (2*Math.PI);
	y = 0.5 - Math.asin(n.e(2))/Math.PI;
	return [x,y]
}

function get_sphere_normal_from_normal_map(texture_map, orig_n, u , v) {
	orig_n = $V2([orig_n.e(1),orig_n.e(2),orig_n.e(3)]);
	var color = get_tga_bilinear_rgb_by_x_y(texture_map, u, v);
	t = 2*color.e(1)-1; b = 2*color.e(2)-1; n =  2*color.e(3)-1;
	normal_tan_plane_coords = $V2([t,b,n]);
	
	// take the 'x' plane vector (compliment to 'y')
	var tan = $V2([0,1,0]).cross(orig_n).toUnitVector();
	// take the compliment
	var cotan = tan.cross(orig_n).toUnitVector();
	
	// transform coords from tangent plane to world plane (inverse the matrix = transpose since it's orthogonal)
	var coord_mat = $M([tan.elements,cotan.elements,orig_n.elements]).transpose();
	res = coord_mat.multiply(normal_tan_plane_coords);
	return $V(res.e(1), res.e(2), res.e(3));
}

function B1_shoot_rays_recursively(pos, ray_d, in_air, recursion_depth) {
	var ret = shoot_ray_and_calc_phong(pos, ray_d); var E  = ret[0]; var color = E;
	if (recursion_depth >= 1) {
		if (E.e(1)!=0 || E.e(2)!=0 || E.e(3)!=0) {
			var coll_obj = ret[1]; var coll_loc = ret[2]; var obj_n = ret[3];
			var reflected_r = obj_n.multiply(2*(ray_d.dot(obj_n))).subtract(ray_d); reflected_r.toUnitVectorN();
			var reflection = B1_shoot_rays_recursively(coll_loc, reflected_r, in_air, recursion_depth-1).multiply(0.5);
			var refraction = $V(0,0,0);
			if (coll_obj.b1_refraction) {
				var idx_r = 1/coll_obj.b1_ref_idx; // idx ratio
				if (~in_air) {
					idx_r = 1/idx_r;
				}
				//refract_ray = -idx_r*(ray - (ray dot norml)*normal) -normal*sqrt(1-idx_r^2(1-(ray dot normal)^2))
				r_dot_n = ray_d.dot(obj_n);
				var refract_ray_l =  ray_d.subtract(obj_n.multiply(r_dot_n)).multiply(-idx_r);
				var refract_ray_r = obj_n.multiply(-1*Math.sqrt(1-Math.pow(idx_r,2)*(1-Math.pow(r_dot_n,2))));
				var refract_ray = refract_ray_l.subtract(refract_ray_r);
				refraction = B1_shoot_rays_recursively(coll_loc, refract_ray, ~in_air, recursion_depth-1).multiply(0.5);
			}
			color = E.multiply(0.5).add(reflection).add(refraction);
		}
		}
	return color;
}

function find_intersections(init_pos, ray_d, objs) {
	if (G.Use_Octree) {
		objs = G.octree.getIntersectionObjects(init_pos, ray_d);
	}
	var int_obj = null; var int_n = null;
	if (objs.length > 0) {
		var i; var intersected_objs =[]; var int_loc = null; var int_idx = null; var int_min_loc_intersec = null;
		// TODO: hack but I don't need there's something prettier needed
		var int_method = intersect_sphere;
		if (objs[0].Type == "Triangle") { 
			int_method = intersect_triangle;
		} else if (objs[0].Type == "Hemisphere") {
			return find_intersections_boolean(init_pos,ray_d, objs);
		}
		
		for (i=0; i < objs.length; i++ ) {
			var obj_intersec = int_method(init_pos, ray_d, objs[i]);
			if (obj_intersec != null && obj_intersec > G.EPSILON) {
				if (int_min_loc_intersec == null || obj_intersec < int_min_loc_intersec) {
					int_min_loc_intersec = obj_intersec;
					int_loc = init_pos.add(ray_d.multiply(obj_intersec));
					int_idx = i;
				}
			}
		}
		
		// yuck!
		if (int_loc!=null) {
			int_obj = objs[int_idx];
			var int_n = get_normal(int_obj, int_loc);
			if (int_n.dot(ray_d) > 0) int_n.multiplyN(-1);
		}
	}
	return [int_obj, int_loc, int_n];
}

function find_intersections_boolean(init_pos,ray_d, objs) {
	int_obj = null; int_loc = null; int_n = null;
	spheres_int = objs[1];
	
	// Bha!
	hem_obj = objs[0];
	// e.g = 0 (but leave a place for numerical error), this means they are parallel, there are 2 cases, but I'm lazy. One will suffice here.
	plane_n = hem_obj.Plane_Normal;
	ray_dot_plane_n = ray_d.dot(hem_obj.Plane_Normal); var hem_t = null;
	if (Math.abs(ray_dot_plane_n) > 0.0001) {
		// Sphere intersection, again.
		var v = init_pos.subtract(hem_obj.Center); var v_dot_d = v.dot(ray_d); var r = hem_obj.Radius;
		var delta = Math.pow(v_dot_d,2)-(v.dot(v) - r*r);
		if (delta > 0) { // no real solution -> no collision
			var sqrt_val = Math.sqrt(delta);

			var t_first = -v_dot_d - sqrt_val; var t_second = -v_dot_d + sqrt_val;
			var plane_t = (hem_obj.Plane_Pt.subtract(init_pos).dot(hem_obj.Plane_Normal))/ray_dot_plane_n;
			if (t_first > 0.001 && t_second > 0.001) {
				if (plane_t > t_first && plane_t < t_second) {
					int_loc = init_pos.add(ray_d.multiply(plane_t)); int_obj = YellowPlane; int_n = hem_obj.Plane_Normal;
					var hem_t = plane_t;
				} else if (t_first > plane_t){
					int_loc = init_pos.add(ray_d.multiply(t_second)); int_obj = LargeRedSphere; int_n = get_normal(int_obj, int_loc);
					var hem_t = t_second;
				}
			}
		}
	}
	
	var sph_int1 = intersect_sphere(init_pos, ray_d, objs[1].Sphere1); var sph_int2 = intersect_sphere(init_pos, ray_d, objs[1].Sphere2);
	if (sph_int1 && sph_int2) {
		var obj_idx = 1; var sph_int_obj = objs[1].Sphere1;
		if (sph_int2 < sph_int1) {
			sph_int1 = sph_int2; sph_int_obj = objs[1].Sphere2;
		}
		if (hem_t == null || sph_int1 < hem_t) {
			int_loc = init_pos.add(ray_d.multiply(sph_int1)); int_n = get_normal(sph_int_obj, int_loc); int_obj = sph_int_obj;
		}
	}
	
	return [int_obj, int_loc, int_n];
}

function intersect_sphere(init_pos, d, sphere) {

	var v = init_pos.subtract(sphere.Center); var v_dot_d = v.dot(d); var r = sphere.Radius;
	var delta = Math.pow(v_dot_d,2)-(v.dot(v) - r*r);
	if (delta < 0) {return null;} // no real solution -> no collision
	var sqrt_val = Math.sqrt(delta);
	
	var t1 = -v_dot_d + sqrt_val; var t2 = -v_dot_d - sqrt_val;
	
	if (t2 > 0.001) {return t2;}
	if (t1 > 0.001) {return t1;}
	return null;
}

function intersect_triangle(init_pos, d, tr) {
	return tr.intersect(init_pos, d);
}
function intersect_ellipsoid(init_pos, ray_d, elps) {
	
	M = $V(1/elps.Radii.e(1),1/elps.Radii.e(2),1/elps.Radii.e(3));
	V1 = scal_vec_mult(ray_d,M);
	PM = scal_vec_mult(init_pos, M);
	P1 = PM.subtract(scal_vec_mult(elps.Center, M));
	
	P1_dot_V1 = P1.dot(V1);
	
	var a = V1.dot(V1); var b = 2*P1_dot_V1; var c = P1.dot(P1) - 1;
	
	var delta = Math.pow(b,2) - 4*a*c;
	if (delta < 0) {return null;} // no real solution -> no collision
	
	var sqrt_val = Math.sqrt(delta);
	var t1 = (-b+sqrt_val)/(2*a); var t2 = (-b-sqrt_val)/(2*a);
	
	if (t2 > 0.001) {return t2;}
	if (t1 > 0.001) {return t1;}
	
	return null;
}

function intersect_elliptic_cylinder(init_pos, ray_d, ell_cyl) {
	M = $V(1/ell_cyl.Radii_x,0,1/ell_cyl.Radii_z);
	
	V1 = scal_vec_mult(ray_d,M);
	P1 = scal_vec_mult(init_pos, M);
	P1_dot_V1 = P1.dot(V1);
	
	var a = V1.dot(V1); var b = 2*P1_dot_V1; var c = P1.dot(P1) - 1;
	
	var delta = Math.pow(b,2) - 4*a*c;
	if (delta < 0) {return null;} // no real solution -> no collision
	
	var sqrt_val = Math.sqrt(delta);
	var t1 = (-b+sqrt_val)/(2*a); var t2 = (-b-sqrt_val)/(2*a);
	
	if (t2 > 0.001) {return t2;}
	if (t1 > 0.001) {return t1;}
	
	return null;
}

function get_normal(int_obj, int_loc) {
	var n = null;
	if (int_obj.Type == "Sphere") {
		n = int_loc.subtract(int_obj.Center);
	} else if (int_obj.Type == "EllipticCylinder") {
		// TODO: hackier sine center is at 0,0
		n = $V((2*int_loc.e(1))/Math.pow(int_obj.Radii_x,2), 0, (2*int_loc.e(3))/Math.pow(int_obj.Radii_z,2));
	} else if (int_obj.Type == "Triangle") {
		n = int_obj.getRayPtNormal(int_loc);
	}
	
	if (n!=null) {
		n.toUnitVectorN();
	}
	return n;
}

function scal_vec_mult(v1,v2) {
	v_x = v1.e(1) * v2.e(1); v_y = v1.e(2) * v2.e(2); v_z = v1.e(3) * v2.e(3);
	return $V(v_x,v_y,v_z);
}

function get_tan_deg(deg) {
   var rad = deg * Math.PI/180;
   return Math.tan(rad);
}

function get_rgb_by_x_y(myTexture, col, row) {
	var id = 3 * (row * myTexture.header.width + col);

	var r = myTexture.image[id + 2] / 255.0;
	var g = myTexture.image[id + 1] / 255.0;
	var b = myTexture.image[id + 0] / 255.0;
	return $V(r,g,b);
}

function get_tga_bilinear_rgb_by_x_y(texture, col, row) {
	// fetch a bilinearly filtered texel
    var u = col * texture.header.width;
    var v = row * texture.header.height;

    // integer pixel values surrounding the exact value
    var u1 = Math.floor(u);
    var v1 = Math.floor(v);
    var u2 = (u1 + 1);
    var v2 = (v1 + 1);
    // calculate fractional parts of u and v
    var frac_u = u - Math.floor (u);
    var frac_v = v - Math.floor (v);
    // calculate weights
    var w1 = (1 - frac_u) * (1 - frac_v);
    var w2 = frac_u * (1 - frac_v);
    var w3 = (1 - frac_u) * frac_v;
    var w4 = frac_u *  frac_v;

    // scale and sum surrounding pixels
    var p1 = get_rgb_by_x_y(texture, u1, v1).multiply(w1);
    var p2 = get_rgb_by_x_y(texture, u2, v1).multiply(w2);
    var p3 = get_rgb_by_x_y(texture, u1, v2).multiply(w3);
    var p4 = get_rgb_by_x_y(texture, u2, v2).multiply(w4);
    return p1.add(p2).add(p3).add(p4);
}

function getSoftShadow(ray_d, intersection, radius, samples_num) {
	// Get a base for the disc (the light turns to 0,0,0)
	var upDirection = G.light_p.cross($V(1,0,0)).toUnitVector();
    var rightDirection = G.light_p.cross(upDirection).toUnitVector();

    var subPositions = new Array();

	// Generate "samples_num" pt light sources
	for (var i = 0; i < samples_num; i++) {
		do {
		var x = Math.random()*2 - 1;
		var y = Math.random()*2 - 1;
		subPositions[i] = G.light_p.add (upDirection.multiply(y)).add (rightDirection.multiply(x));
		} while (subPositions[i].subtract(G.light_p).modulus() > radius);
	}

    var shadow_light_scalar = 1;

    for (var i = 0; i < samples_num; i++) {
        var monte_shadow_ray = subPositions[i].subtract(intersection).toUnitVector();

		var light_intersection = find_intersections(intersection, monte_shadow_ray, G.SceneObjects); var intersec_obj = light_intersection[0];
        if (intersec_obj !== null)  {
			shadow_light_scalar -= 1/samples_num;
		}
    }

    return shadow_light_scalar;
}