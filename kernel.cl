#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//To render the random_scene, choose the background from the commented options,
//Then change the world and mats array length to 500.
__constant int MAX_DEPTH = 20;
__constant int NUM_POINTS = 100;
__constant int USE_HAMM = 0;
__constant int NUM_U = 8;
__constant int NUM_V = 16;
__constant int NUM_ACTIONS = 128;
__constant int RL_ON=1;  //0 for without RL and 1 for rendering with RL ON
__constant double FUZZ = 0.01;
__constant double3 background = (double3){0.0,0.0,0.0};
// __constant double3 background = (double3){0.05,0.05,0.05};
// __constant double3 background = (double3){0.68, 0.84, 0.9};

// Define C++ Classes as OpenCL structs
typedef struct _cl_tag_sObject {
	double x0,x1,y0,y1,z0,z1,k,radius;
	int m_type;
	int off;
	int numPoints;
} sObject;
typedef struct _cl_tag_Ray {
	double3 a;
	double3 b;
} Ray;
typedef struct _cl_tag_Camera {
	double3 lookFrom;
	double3 lookAt;
	double3 viewUp;
	double aperture;
	double Fov;
	double focus_dist;
} Camera;
typedef struct _cl_tag_Material {
	double3 m_vColor;
	bool flip_face;
	int m_MType;
} Material;

typedef struct _cl_tag_HitRecord {
	double3 p;
	double3 normal;
	Material m_curmat;
	double t,u,v;
	bool front_face;
	int objId;
} HitRecord;
void set_face_normal(Ray r,double3 outward_normal,HitRecord *rec){
	rec->front_face = dot(r.b,outward_normal) < 0;
	if (rec->front_face)
		rec->normal = outward_normal;
	else
		rec->normal = -outward_normal;
}
// Define math functions
double3 UnitVector(double3 v) {
	return v / length(v);
}
double3 PointAtParameter(Ray r, double t) { 
	return r.a + t * r.b; 
}
double3 InvDir(const Ray r) {
	return 1 / r.b;
}
double SquaredLength(double3 m_dE) {
	return m_dE.x * m_dE.x + m_dE.y * m_dE.y + m_dE.z * m_dE.z;
}
double SquaredLength2(double2 m_dE) {
	return m_dE.x * m_dE.x + m_dE.y * m_dE.y;
}
uint MWC64X(uint2 *state)
{
    enum { A=4294883355U};
    uint x=(*state).x, c=(*state).y;  // Unpack the state
    uint res=x^c;                     // Calculate the result
    uint hi=mul_hi(x,A);              // Step the RNG
    x=x*A+c;
    c=hi+(x<c);
    *state=(uint2){x,c};            // Pack the state back up
    return res;                       // Return the next result
}
double myrand(uint2* seed_ptr) //Random Double
{
	uint MAX_INT=0;
	MAX_INT--;
	uint rInt = MWC64X(seed_ptr);
    double ans = ((double)rInt)/MAX_INT;
    return ans;
 }
 double3 random_in_unit_disk(uint2* seed_ptr) {
    while (true) {
        double3 p = (double3){2*myrand(seed_ptr)-1,2*myrand(seed_ptr)-1,0};
        if (SquaredLength(p) >= 1) continue;
        return p;
    }
}
double2 get_sphere_uv(double3 p) {    //texture map of sphere
	double2 ans;
	double phi = atan2(p.z, p.x);
	double theta = asin(p.y);
	ans.x = 1-(phi + M_PI) / (2*M_PI);
	ans.y = (theta + M_PI/2) / M_PI;
	return ans;
}
double2 get_obj_uv(sObject obj,double3 p){  //texture map of all objects
	double2 uv;
	if (obj.m_type == 1){
		uv.x = (p.x-obj.x0)/(obj.x1-obj.x0);
		uv.y = (p.y-obj.y0)/(obj.y1-obj.y0);
	}
	else if (obj.m_type == 2){
		uv.x = (p.x-obj.x0)/(obj.x1-obj.x0);
		uv.y = (p.z-obj.z0)/(obj.z1-obj.z0);
	}
	else if (obj.m_type == 3){
		uv.x = (p.y-obj.y0)/(obj.y1-obj.y0);
		uv.y = (p.z-obj.z0)/(obj.z1-obj.z0);
	}
	else if (obj.m_type == 4){
		double3 center = (double3){obj.x0,obj.y0,obj.z0};
		uv = get_sphere_uv((p-center)/obj.radius);
	}
	return uv;
}
int getClosestIndexHamm(global const double2* points, sObject obj,double3 p){  //Linear search for closest point
	double2 uv = get_obj_uv(obj,p);
    int ans=-1;
    double mindis=0;
    for (int i=0;i<obj.numPoints;i++){
		double2 temp = points[i+obj.off] - uv;
        double dis = sqrt(SquaredLength2(temp));
        if (ans==-1){
            ans = i;
            mindis = dis;
        }
        else if (dis < mindis){
            mindis = dis;
            ans = i;
        }
    }
    return ans;
}
int getClosestIndexGrid(sObject obj,double3 p){  //Linear search for closest point
	double2 uv = get_obj_uv(obj,p);
    int xid = min((int)(uv.y*obj.numPoints),obj.numPoints-1);
	int yid = min((int)(uv.x*obj.numPoints),obj.numPoints-1);
	int ans = yid*obj.numPoints + xid;
    return ans;
}
double3 random_cosine_direction(uint2* seed_ptr){
	double r1 = myrand(seed_ptr);
	double r2 = myrand(seed_ptr);
	double z = sqrt(1-r2);
	double phi = 2*M_PI*r1;
	double x = cos(phi)*sqrt(r2);
	double y = sin(phi)*sqrt(r2);
	return (double3){x,y,z};
}

// Define external functions
double myabs(double x){
	if (x<0)
		x=-x;
	return x;
}
double3 Reflect(const double3 v, const double3 n) {
	return (double3)(v - 2 * dot(v, n)*n);
}
double3 Refract(const double3 uv, const double3 n, double etai_over_etat) {
	double cos_theta = min(dot(-uv, n), 1.0);
    double3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    double3 r_out_parallel = -sqrt(myabs(1.0 - SquaredLength(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}
double schlick(double cosine, double ref_idx) {
	double d0 = (1 - ref_idx) / (1 + ref_idx);
	d0 = d0 * d0;
	return d0 + (1 - d0)*pow((1 - cosine), 5);
}
double3 get_rand_dir(int idx_Q,uint2* seed_ptr, bool center){  //Random direction according to idxQ if center is false
	int numu = NUM_U, numv = NUM_V;
	double delu = 1.0/numu, delv = 1.0/numv;
	int uid = idx_Q/numv, vid = idx_Q%numv;
	double v_min = vid*delv, u_min = uid*delu;
	double r1,r2;
	if (center){
		r1 = 0.5;
		r2 = 0.5;
	}
	else{
		r1 = myrand(seed_ptr);
		r2 = myrand(seed_ptr);
	}
	double u = u_min + delu*r1;
	double v = v_min + delv*r2;
	double phi = 2*M_PI*v, temp = sqrt(1-u*u);
	double3 rand_dir = (double3){temp*cos(phi),temp*sin(phi),u};
	return rand_dir;
}
double3 get_direction(double3 normal,int idx_Q,uint2* seed_ptr, bool center){
	double3 axis[3];
	axis[2] = UnitVector(normal);
	double3 a;
	if (myabs(axis[2].x)>0.9)
		a=(double3){0,1,0};
	else
		a=(double3){1,0,0};
	axis[1] = UnitVector(cross(axis[2],a));
	axis[0] = cross(axis[2],axis[1]);

	double3 rand_dir;
	if (RL_ON==1)
		rand_dir = get_rand_dir(idx_Q,seed_ptr,center);
	else
		rand_dir = random_cosine_direction(seed_ptr);
	
	double3 direction = rand_dir.x*axis[0] + rand_dir.y*axis[1] + rand_dir.z*axis[2];
	return UnitVector(direction);
}
bool scatter(Material *mat, const Ray r_in, const HitRecord *rec, double3 *attenuation, double *pdf_val, Ray *scattered, uint2* seed_ptr, int idx_Q) { //combined scatter function
	int gid = get_global_id(0);
	// Lambertian
	if (mat->m_MType == 0) {
		double3 direction = get_direction(rec->normal,idx_Q,seed_ptr,false);
		*scattered = (Ray){rec->p,direction};
		*attenuation = mat->m_vColor;
		*pdf_val = dot(UnitVector(rec->normal),scattered->b)/M_PI;
		return true;
	}
	// Light
	else if (mat->m_MType == 1) {
		return false;
	}
	// Metal
	else if (mat->m_MType == 2) {
		double3 reflected = Reflect(UnitVector(r_in.b),rec->normal);

		double3 rand_dir;
		if (RL_ON==1)
			rand_dir = get_rand_dir(idx_Q,seed_ptr,false);
		else
			rand_dir = random_cosine_direction(seed_ptr);

		double fuzz = FUZZ;
		*scattered = (Ray) {rec->p, UnitVector(reflected+fuzz*rand_dir)};
		*attenuation = mat->m_vColor;
		*pdf_val = 1;
		return true;
	}
	//Dielectric
	else if (mat->m_MType == 3) {
		*attenuation = mat->m_vColor;
		double ref_idx = 1.5;
		double etai_over_etat = rec->front_face ? (1.0 / ref_idx) : ref_idx;
		double3 unit_direction = UnitVector(r_in.b);
		double cos_theta = min(dot(-unit_direction, rec->normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

		double3 rand_dir;
		if (RL_ON==1)
			rand_dir = get_rand_dir(idx_Q,seed_ptr,false);
		else
			rand_dir = random_cosine_direction(seed_ptr);
		double fuzz = 0.01;

		if (etai_over_etat * sin_theta > 1.0 ) {
			double3 reflected = Reflect(unit_direction, rec->normal);
			*scattered = (Ray) {rec->p, UnitVector(reflected+fuzz*rand_dir)};
			*pdf_val = 1;
			return true;
		}
		double reflect_prob = schlick(cos_theta, etai_over_etat);
		if (myrand(seed_ptr) < reflect_prob)
		{
			double3 reflected = Reflect(unit_direction, rec->normal);
			*scattered = (Ray) {rec->p, UnitVector(reflected+fuzz*rand_dir)};
			*pdf_val = 1;
			return true;
		}
		double3 refracted = Refract(unit_direction, rec->normal, etai_over_etat);
		*scattered = (Ray) {rec->p, refracted};
		*pdf_val = 1;
		return true;
	}
	printf("Hi1\n");
	return false;
}
bool hit(const sObject obj, const Material m, const Ray r, double t0, double t1, HitRecord *rec) {  //combined hit function
	if (obj.m_type==1) {
		if (obj.radius!=0){
			sObject tobj = obj;
			double objangle = obj.radius;
			double angle = objangle*M_PI/180.0;
			double cos_theta = cos(angle), sin_theta = sin(angle);
			tobj.radius = 0;
			double3 origin = r.a;
			double3 direction = r.b;
			// origin.x = origin.x - obj.x1;
			// origin.z = origin.z - obj.z1;
			origin.x = cos_theta*origin.x - sin_theta*origin.z;
			origin.z = sin_theta*origin.x + cos_theta*origin.z;
			direction.x = cos_theta*direction.x - sin_theta*direction.z;
			direction.z = sin_theta*direction.x + cos_theta*direction.z;
			// origin.x = origin.x + obj.x1;
			// origin.z = origin.z + obj.z1;
			Ray rotated_r = (Ray) {origin, direction};
			if (!hit(tobj,m,rotated_r,t0,t1,rec)){
				return false;
			}
			double3 p = rec->p;
			double3 normal = rec->normal;
			p.x = cos_theta*p.x + sin_theta*p.z;
			p.z = -sin_theta*p.x + cos_theta*p.z;
			normal.x = cos_theta*normal.x + sin_theta*normal.z;
			normal.z = -sin_theta*normal.x + cos_theta*normal.z;
			rec->p = p;
			rec->normal = normal;
			return true;
		}
		double t = (obj.k - r.a.z)/r.b.z;
		if (t<t0 || t>t1)
			return false;
		double3 p_at = PointAtParameter(r,t);
		double x = p_at.x;
		double y = p_at.y;
		if (x < obj.x0 || x > obj.x1 || y < obj.y0 || y > obj.y1)
        	return false;
		double2 uv = get_obj_uv(obj,p_at);
		rec->u = uv.x;
		rec->v = uv.y;
		rec->t = t;
		set_face_normal(r,(double3){0,0,1},rec);
		rec->m_curmat = m;
		rec->p = p_at;
		return true;
	}
	else if (obj.m_type==2){
		double t = (obj.k - r.a.y)/r.b.y;
		if (t<t0 || t>t1)
			return false;
		double3 p_at = PointAtParameter(r,t);
		double x = p_at.x;
		double z = p_at.z;
		if (x < obj.x0 || x > obj.x1 || z < obj.z0 || z > obj.z1)
        	return false;
		double2 uv = get_obj_uv(obj,p_at);
		rec->u = uv.x;
		rec->v = uv.y;
		rec->t = t;
		set_face_normal(r,(double3){0,1,0},rec);
		rec->m_curmat = m;
		rec->p = p_at;
		return true;
	}
	else if (obj.m_type==3){
		double t = (obj.k - r.a.x)/r.b.x;
		if (t<t0 || t>t1)
			return false;
		double3 p_at = PointAtParameter(r,t);
		double y = p_at.y;
		double z = p_at.z;
		if (y < obj.y0 || y > obj.y1 || z < obj.z0 || z > obj.z1)
        	return false;
		double2 uv = get_obj_uv(obj,p_at);
		rec->u = uv.x;
		rec->v = uv.y;
		rec->t = t;
		set_face_normal(r,(double3){1,0,0},rec);
		rec->m_curmat = m;
		rec->p = p_at;
		return true;
	}
	else if (obj.m_type==5){
		double3 center = (double3){obj.x0,obj.y0,obj.z0};
		double3 oc = r.a-center;
		double a = SquaredLength(r.b);
		double half_b = dot(oc,r.b);
		double c = SquaredLength(oc)-obj.radius*obj.radius;
		double discriminant = half_b*half_b - a*c;
		if (discriminant>0){
			double root = sqrt(discriminant);
			double temp = (-half_b-root)/a;
			bool poss=false;
			if (temp<t1 && temp>t0)
				poss=true;
			else{
				temp = (-half_b + root) / a;
				poss=temp < t1 && temp > t0;
			}
			if (poss){
				rec->t = temp;
				rec->p = PointAtParameter(r,temp);
				double3 outward_normal=(rec->p-center)/obj.radius;
				set_face_normal(r,outward_normal,rec);
				double2 uv = get_sphere_uv(outward_normal);
				rec->u = uv.x;
				rec->v = uv.y;
				rec->m_curmat = m;
				return true;
			}
			return false;
		}
		return false;
	}
	printf("Hi\n");
	return false;
}
bool worldHit(const sObject *x, const Material *m, int ObjLen, const Ray r, HitRecord *rec) {
	HitRecord temp_rec;
	bool hitAnything = false;
	double closestSoFar = DBL_MAX;
	for (int i = 0; i < ObjLen; i++) {
		if (hit(x[i], m[i], r, 0.001, closestSoFar, &temp_rec)) {
			hitAnything = true;
			closestSoFar = temp_rec.t;
			*rec = temp_rec;
			rec->objId = i;
		}
	}
	return hitAnything;
}
double3 getEmitted(Material *mat,Ray r,HitRecord rec){  //emitted only for light material
	if (mat->m_MType!=1){
		return (double3){0,0,0};
	}
	else{
		bool cond = rec.front_face;
		if (mat->flip_face)
			cond = !rec.front_face;
		if (cond){
			return mat->m_vColor;
		}
		else{
			return (double3){0,0,0};
		}
	}
	printf("Hi\n");
	return (double3){0,0,0};
}
double double3_max(double3 val){  //max over R,G,B
	double ans = val.x;
	if (val.y>ans) 
		ans=val.y;
	if (val.z>ans) 
		ans=val.z;
	return ans;
}
double scattering_pdf(Material *mat, Ray r_in, HitRecord rec,Ray scattered){  //cosine pdf for scattering
	if (mat->m_MType == 0) {
		double cosine = dot(rec.normal,UnitVector(scattered.b));
		if (cosine<0)
			return 0;
		else
			return cosine/M_PI;
	}
	else if (mat->m_MType==1){
		return 0;
	}
	printf("Hi\n");
	return 0;
}
double3 Color(const Ray *start,sObject *world, Material* mats, const int ObjLen, int max_depth, int* count, uint2* seed_ptr, global const double2* hammerselyPoints,global double* qtable) {
	double3 ans = (double3){1,1,1}, lastAttenuation;
	Ray r = *start;
	int idxQ = 0,prevOffset, nA = NUM_ACTIONS;
	for (int depth=0;depth<=max_depth;depth++){
		if (depth >= max_depth){
			ans = ans*(double3){0,0,0};
			break;
		}
		HitRecord rec;
		if (!worldHit(world,mats,ObjLen,r,&rec)){
			ans = ans * background;
			count[0]+=depth;
			break;
		}
		Ray scattered;
		double3 emitted = getEmitted(&rec.m_curmat,r,rec);
		int idxCurr=0,currOff=0;
		double Sum=0;
		if (RL_ON==1){
			if (USE_HAMM==1)
				idxCurr = getClosestIndexHamm(hammerselyPoints,world[rec.objId],rec.p);
			else
				idxCurr = getClosestIndexGrid(world[rec.objId],rec.p);
			currOff = (world[rec.objId].off + idxCurr)*nA;
			if (depth>0){
				double3 update;
				if (rec.m_curmat.m_MType==1)
					update = emitted;
				else{
					int qmaxid = currOff;
					for (int i=currOff;i<currOff+nA;i++){
						if (qtable[i]>qtable[qmaxid])
							qmaxid = i;
					}
					update = lastAttenuation*qtable[qmaxid];
				}
				//This region is for experimenting with Expected Sarsa
				// else {
				// 	int offset = currOff;
				// 	double sval=0,sumval=0;
				// 	for (int i=offset;i<offset+nA;i++) {
				// 		double3 diri = get_direction(rec.normal,i-offset,seed_ptr,true);
				// 		double brdf = dot(rec.normal,diri);
				// 		sval = sval + brdf*qtable[i];
				// 		sumval = sumval + brdf;
				// 	}
				// 	update = (sval/sumval)*lastAttenuation;
				// }
				double lr = 0.2;
				double update_val = double3_max(update);
				qtable[prevOffset] = (1-lr)*qtable[prevOffset] + lr*update_val;
			}
			double scatter_cdf[128];
			for (int i=0;i<nA;i++){
				Sum=Sum+qtable[i+currOff];
				scatter_cdf[i]=Sum;
			}
			double temprand = Sum*myrand(seed_ptr);
			idxQ = 0;
			while(idxQ<nA){
				if (temprand<=scatter_cdf[idxQ] || idxQ==nA)
					break;
				idxQ = idxQ + 1;
			}
		}
		
		double3 albedo;
		double pdf_val;
		if (!scatter(&rec.m_curmat,r,&rec,&albedo,&pdf_val,&scattered, seed_ptr, idxQ)){
			ans = ans * emitted;
			break;
		}
		if (RL_ON==1)
			pdf_val = (qtable[currOff+idxQ]/Sum)*nA/(2*M_PI);
		double3 attenuation;
		if (rec.m_curmat.m_MType==2 || rec.m_curmat.m_MType==3){
			attenuation = albedo;
			pdf_val = 1;
		}
		else{
			double spdf = scattering_pdf(&rec.m_curmat,r,rec,scattered);
			attenuation = albedo * spdf;
		}
		ans = ans * attenuation / pdf_val;
		r = scattered;
		if (RL_ON==1){
			lastAttenuation = attenuation;
			prevOffset = currOff + idxQ;
		}
	}
	return ans;
}
Ray getRay(double s, double t, int2 dims, Camera cam,uint2* seed_ptr) { //get Ray from pixel according to (s,t) in pixel

	double dHalfHeight = tan(cam.Fov*M_PI / 360);
	double dHalfWidth = ((double)dims.x / dims.y) * dHalfHeight;
	double dFocusDist = length(cam.lookFrom - cam.lookAt);
	// double dFocusDist = cam.focus_dist;
	double3 vW = UnitVector(cam.lookFrom - cam.lookAt);
	double3 vU = UnitVector(cross(cam.viewUp, vW));
	double3 vV = cross(vW, vU);
	double3 vOrigin = cam.lookFrom;
	double3 vLowerLeftCorner = vOrigin - (dHalfWidth * dFocusDist * vU) - (dHalfHeight * dFocusDist * vV) - (dFocusDist * vW);
	double3 vHorizontal = 2 * dHalfWidth*dFocusDist*vU;
	double3 vVertical = 2 * dHalfHeight*dFocusDist*vV;
	double3 vRD = (cam.aperture / 2) * random_in_unit_disk(seed_ptr);
	double3 vOffset = vU * vRD.x + vV * vRD.y;
	return (Ray) { (double3)(vOrigin + vOffset), (double3)(vLowerLeftCorner + (s * vHorizontal) + (t * vVertical) - vOrigin - vOffset) };
}

kernel void Render(global double4 *pixel, global int2 *dims, global const double16 *cam, global const double8 *list, global const int *listLen, global const double4 *materials, global const double2* hammerselyPoints,global const int* ssp, global const int* hamOff,global double* qtable) {

	int gid = get_global_id(0); // Current ray in image
	int raysPerPixel = ssp[0];
	int gsize = get_global_size(0); // Number of rays in pixel
	int2 dim = dims[0]; // Image Dimensions
	int ObjLen = listLen[0]; // # objects in list
	int wgNum = get_group_id(0);
	
	int j = floor((double)(1 + (gid / dim.x))); // Current Y
	int i = gid - ((j - 1) * dim.x); // Current X
	// Object list initialized from kernel struct, max objects in image defined by array size
	sObject world[20];
	Material mats[20];
	Camera camx[1];
	int tot_points=0;
	for (int i = 0; i < ObjLen; i++) {
		world[i].m_type = (int)list[i].s7;
		if (USE_HAMM==1){
			if (i==0){
				world[i].off = 0;
				world[i].numPoints = hamOff[i];
			}
			else{
				world[i].off = hamOff[i-1];
				world[i].numPoints = hamOff[i] - hamOff[i-1];
			}
		}
		else{
			world[i].off = tot_points;
			world[i].numPoints = hamOff[i];
			tot_points = tot_points + hamOff[i]*hamOff[i];
		}
		if (world[i].m_type == 1){
			world[i].x0 = list[i].s0;
			world[i].x1 = list[i].s1;
			world[i].y0 = list[i].s2;
			world[i].y1 = list[i].s3;
			world[i].z0 = world[i].z1 = world[i].k = list[i].s4;
			world[i].radius = list[i].s6;  //rotation amount stored in radius for rectangles
		}
		else if (world[i].m_type == 2){
			world[i].x0 = list[i].s0;
			world[i].x1 = list[i].s1;
			world[i].z0 = list[i].s2;
			world[i].z1 = list[i].s3;
			world[i].y0 = world[i].y1 = world[i].k = list[i].s4;
			world[i].radius = list[i].s6;
		}
		else if (world[i].m_type == 3){
			world[i].y0 = list[i].s0;
			world[i].y1 = list[i].s1;
			world[i].z0 = list[i].s2;
			world[i].z1 = list[i].s3;
			world[i].x0 = world[i].x1 = world[i].k = list[i].s4;
			world[i].radius = list[i].s6;
		}
		else if (world[i].m_type == 5){
			world[i].x0 = list[i].s0;
			world[i].y0 = list[i].s1;
			world[i].z0 = list[i].s2;
			world[i].radius = list[i].s3;
		}

		mats[i].m_vColor = (double3){ materials[i].s0, materials[i].s1, materials[i].s2 };
		mats[i].m_MType = (int)(materials[i].s3);
		mats[i].flip_face = mats[i].m_MType == 1 && world[i].m_type != 5;
	}
	camx[0].lookFrom = (double3) { cam[0].s0, cam[0].s1, cam[0].s2 };
	camx[0].lookAt = (double3) { cam[0].s3, cam[0].s4, cam[0].s5 };
	camx[0].viewUp = (double3) { cam[0].s6, cam[0].s7, cam[0].s8 };
	camx[0].Fov = cam[0].s9;
	camx[0].aperture = cam[0].sa;
	camx[0].focus_dist = cam[0].sb;
	
	double3 col = (double3)(0);

	int count[1];
	count[0]=0;
	int max_depth = MAX_DEPTH;
	uint2 SEED;
	SEED.x = gid;
	SEED.y = SEED.x+1;
	int iters = 1;
	if (RL_ON==1 && gid<dim.x*dim.y/10)
		iters = 2 + NUM_POINTS/50;
		// iters = 1;
	for (int iter=0;iter<iters;iter++){
		col = (double3)(0);
		for (int s = 0; s < raysPerPixel; s++) {
			double u = (double)(i + myrand(&SEED)) / (double)dim.x;
			double v = (double)(j + myrand(&SEED)) / (double)dim.y;
			Ray r = getRay(u, v, dim, camx[0],&SEED);
			col += Color(&r, &world, &mats, ObjLen, max_depth, &count, &SEED, hammerselyPoints,qtable);
		}
	}
	// col = sqrt(col / raysPerPixel);             //No need for gamma correction as it is done while writing to image file by CPU
	
	pixel[gid] = (double4)(col, gid);
	double numElements = dim.x*dim.y;
	printf("\r%.2f percent",100*gid/numElements);  //Print progress
}