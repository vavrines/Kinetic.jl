# ============================================================
# Coffee Creamer of Quadrature
# ============================================================

function getdist(v1::Array{Float64,1},v2::Array{Float64,1},geometry="planar")
	if geometry=="planar"
		return norm(v1-v2)
	elseif geometry=="sphere"
		if (v1[1]-v2[1])*(v1[1]-v2[1])+
			(v1[2]-v2[2])*(v1[2]-v2[2])+
			(v1[3]-v2[3])*(v1[3]-v2[3])<1e-5
			#if norm(v1-v2)<1e-15
			return norm(v1-v2)
		end
		#https://en.wikipedia.org/wiki/Great-circle_distance
		#return atan(norm(cross(v1,v2),2)/dot(v1,v2))
		return acos(max(-1.0,v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]))
	else 
		error("Type has to be planar or sphere.")
	end
end

function getangle(B::Array{Float64,1},A::Array{Float64,1},C::Array{Float64,1}, geometry="planar")
	if geometry=="planar"
		u,v = A-B, C-A;
		return acos(dot(u,v)/norm(u,2)/norm(v,2));
	elseif geometry=="sphere"
		c= getdist(B,A,"sphere");
		b =getdist(A,C,"sphere");
		a= getdist(C,B,"sphere");
		if c<1e-10 || b<1e-10 || a<1e-10
			return 0;
		end
		#https://en.wikipedia.org/wiki/Spherical_trigonometry#Cosine_rules_and_sine_rules
		tmp = (cos(a)-cos(b)*cos(c))/(sin(b)*sin(c));
		# It might happen, that due to computations, tmp is not inside [-1,1]
		# but inside [-1-eps,1+eps]
		tmp = max(-1.0,min(1.0,tmp))

		angle = acos(tmp);
		return angle;
	else 
		error("Type has to be planar or sphere.")
	end
end

function getarea(A::Array{Float64,1},B::Array{Float64,1},C::Array{Float64,1}, geometry="planar")
	if geometry=="planar"
		alpha = getangle(B,A,C)
		lb = norm(B-A)
		la = norm(C-A)
		return 0.5*sin(alpha)*lb*la
	elseif geometry=="sphere"
		alpha = getangle(B,A,C,"sphere")
		beta =  getangle(C,B,A,"sphere")
		gamma = getangle(A,C,B,"sphere")
		#https://en.wikipedia.org/wiki/Spherical_trigonometry#Area_and_spherical_excess
		#Excess = max(0.0,alpha+beta+gamma-pi)
		#R = 1; # We always consider the uni sphere
		#return Excess*R*R;
		return alpha+beta+gamma-pi
	else 
		error("Type has to be planar or sphere.")
	end
end

function muphi2xyz!(muphi::Array{Float64,2},xyz::Array{Float64,2})
	n = size(xyz)[1]
	for i=1:n
		xyz[i,1] = sqrt(1-muphi[i,1]^2)*cos(muphi[i,2])
		xyz[i,2] = sqrt(1-muphi[i,1]^2)*sin(muphi[i,2])
		xyz[i,3] =        muphi[i,1]
	end
end

function xyz2muphi!(xyz::Array{Float64,2},muphi::Array{Float64,2})
	n = size(xyz)[1]
	for i=1:n
		muphi[i,1] = xyz[i,3]
		muphi[i,2] =  atan(xyz[i,2],xyz[i,1])    
	end
end

#Check if the projection of p onto the plane defined by p1 p2 p3 is inside the triangle p1 p2 p3
#https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
function isprojectedpointinsidetriangle(p::Array{Float64,1},p1::Array{Float64,1},
	p2::Array{Float64,1} ,p3::Array{Float64,1})
	n1 = cross(p1,p2);
	if n1'p3<0
		if p'n1>0
			return false
		end
	elseif n1'p3>0
		if p'n1<0
			return false
		end
	end

	n2 = cross(p2,p3);
	if n2'p1<0
		if p'n2>0
			return false
		end
	elseif n2'p1>0
		if p'n2<0
			return false
		end
	end

	n3 = cross(p3,p1);
	if n3'p2<0
		if p'n3>0
			return false
		end
	elseif n3'p2>0
		if p'n3<0
			return false
		end
	end
	return true
end

#Compute the rotation that rotates theta around
#the vector u 
function getrotationmatrix!(Rot::Array{Float64,2},u::Array{Float64,1},theta::Float64,order="forward")
	#https://en.wikipedia.org/wiki/Rotation_matrix
	ux = u[1];
	uy = u[2];
	uz = u[3];
	costheta = cos(theta);
	sintheta = sin(theta);

	Rot[1,1] = ux*ux*(1-costheta)+   costheta;
	Rot[1,2] = ux*uy*(1-costheta)-uz*sintheta;
	Rot[1,3] = ux*uz*(1-costheta)+uy*sintheta;

	Rot[2,1] = uy*ux*(1-costheta)+uz*sintheta;
	Rot[2,2] = uy*uy*(1-costheta)+   costheta;
	Rot[2,3] = uy*uz*(1-costheta)-ux*sintheta;

	Rot[3,1] = uz*ux*(1-costheta)-uy*sintheta;
	Rot[3,2] = uz*uy*(1-costheta)+ux*sintheta;
	Rot[3,3] = uz*uz*(1-costheta)+   costheta;
	#niceprint(Rot)
	#if (abs(norm(Rot)-1)>1e-6)
	#	error("No rotation matrix.")
	#end
end

#
#Compute the rotation that rotas alpha around the x axis, beta around the y axis
#and gamma around the z axis.
#
function getrotationmatrix!(Rot::Array{Float64,2},alpha::Float64,beta::Float64,gamma::Float64,order="forward")
	RotX = zeros(3,3)
	RotY = zeros(3,3)
	RotZ = zeros(3,3)
	# set nonzero entries

	RotX[2,2] =  cos(alpha);
	RotX[2,3] = -sin(alpha);
	RotX[3,2] =  sin(alpha);
	RotX[3,3] =  cos(alpha);

	RotY[1,1] =  cos(beta);
	RotY[1,3] =  sin(beta);
	RotY[3,1] = -sin(beta);
	RotY[3,3] =  cos(beta);

	RotZ[1,1] =  cos(gamma);
	RotZ[1,2] = -sin(gamma);
	RotZ[2,1] =  sin(gamma);
	RotZ[2,2] =  cos(gamma);

	if order=="forward"
		Rot = deepcopy(RotZ*RotY*RotX);
	elseif order=="reverse"
		Rot = deepcopy(RotX*RotY*RotY);
	else
		error("Order has to be forward or reverse!")
	end
end

function computeinterpolationweights(v::Array{Float64,1},va::Array{Float64,1},vb::Array{Float64,1},
									 vc::Array{Float64,1})
	#https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf

	weights = zeros(3)

	areac = getarea(v,va,vb,"sphere");
	areab = getarea(v,va,vc,"sphere");
	areaa = getarea(v,vb,vc,"sphere");

	area = areaa+areab+areac
	weights[1] = areaa/area;
	weights[2] = areab/area;
	weights[3] = areac/area;

	return weights

end

function interpolate(v::Array{Float64,1},va::Array{Float64,1},vb::Array{Float64,1},
					 vc::Array{Float64,1},phia::Float64,phib::Float64,phic::Float64)
	w = computeinterpolationweights(v,va,vb,vc)
	return phia*w[1]+phib*w[2]+phic*w[3]
end

function slerp(pt1::Array{Float64,1},pt2::Array{Float64,1},n::Int64,unevenfix=false::Bool)
	if norm(pt1-pt2)<1e-10 # same points
		return repeat(pt1,1,n) # return n copies of that point
	end

	omega =acos(dot(pt1,pt2));
	t = collect(range(0,stop=1,length=n))


	return (sin.((1 .-t)*omega)/sin.(omega))'.*pt1+
	(sin.(t*omega)/sin.(omega))'.*pt2
end

function unique(Points::Array{Float64,2},Triangles::Array{Int64,2})
	# Given points P and a triangulation T, this function removes
	# all duplicate (non unique) entries from P and adjusts the entries of T accordingly.
	# The final points and triangulation are stored in
	# XYZ and Triangulation (i.e.) the class attributes

	nPoints = size(Points)[2]
	nTriangles = size(Triangles)[2]

	map    = collect(1:nPoints)
	unique = ones(Int64,nPoints) 


	#Okay so what is map?
	#map is a vector, where each entry contains the mapping to the ID in the unique vector we want to have.
	#We start with map(i) = i and then iterate through all the points
	#if it turns out, that point 10 == point 4, we will set map(10) = 4 (and map(4)=4)
	#if then also point 20 == point 10 == point 4, we will set map(10) == 4 as well.
	#unique simply stores a boolean as for whether or not the value is unique or a duplicate.
	#This does however keep the first occurrence as unique.
	#Meaning for a vector
	#[1 2 3 3 4 5 3]
	#unique would be
	#[1 1 1 0 1 1 0]
	#
	for i=1:nPoints
		for j=i+1:nPoints
			# apparently calling norm(p1-p2) allocates a lot of
			dist = sqrt(
						(Points[1,i]-Points[1,j])^2+
						(Points[2,i]-Points[2,j])^2+
						(Points[3,i]-Points[3,j])^2)
			if (dist<1e-6) # Points are equal
				map[j] = min(map[j],i); # Take the minimal ID
				unique[j] = 0; # This point is no longer unique
			end
		end
	end

	#We are still not done
	#Consider again the vector from above
	#[1 2 3 3 4 5 3]
	#map was at first
	#[0 1 2 3 4 5 6]
	#and is now
	#[0 1 2 2 4 5 2]
	#meaning that we have to renumber certain entries because
	#we want to get
	#[0 1 2 2 3 4 2] , i.e. the "real count
	#
	uniqueCounter = 1 # The "real" count
	for i=1:nPoints
		if unique[i]==1
			map[i] = uniqueCounter;
			uniqueCounter +=1;
		else
			map[i] = map[map[i]]
		end
	end

	# Now we go through the triangles and update the IDs
	triangulation = zeros(Int64,3,nTriangles) 
	for i=1:nTriangles
		for j=1:3
			idx = Triangles[j,i];
			triangulation[j,i] = min(idx,map[idx])
		end
	end

	#Count unique points
	nQuadPoints=0;
	for i=1:nPoints
		nQuadPoints += (unique[i]==1);
	end

	#Write unique points to XYZ
	xyz = zeros(Float64,3,nQuadPoints);
	counter = 1;
	for i=1:nPoints
		if unique[i]==1
			xyz[:,counter] = Points[:,i]
			counter +=1;
		end
	end
	return xyz,triangulation
end
