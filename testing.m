clc 
clear 

r = [ 42488261.00028026, -50635537.69016445, 0];
v = [-754.37804831, 696.12410414, 0];
mu = 4904869500000.0;

vector_to_elements(r, v, mu)


function [a e i lan w ta] = vector_to_elements(rvec, vvec, mu)
	r = norm(rvec);
	v = norm(vvec); 
	vr = dot(rvec, vvec)/r;
	
	hvec = cross(rvec, vvec); 
	h = norm(hvec);
	
	% Inclination
	i = acos(hvec(3)/h);
	
	% LAN 
	nvec = cross([0 0 1], hvec);
	n = norm(nvec);
	
	if n == 0
		lan = 0;
	else
		lan = acos(nvec(1)/n);
		
		if nvec(2) < 0
			lan = 2*pi - lan;
		end 
	end
	
	% eccentricity 
	evec = (1/mu)*((v^2 - mu/r)*rvec - r*vr*vvec);
	e = norm(evec);
	
	evec/e
	
	% argument of periapsis 
	if n == 0 
		w = acos(dot(evec, [1 0 0])/e)
		
		if evec(3) <= 0
			w = 2*pi - w
		end
	else
		pass 
	end
	
	rad2deg(w)
end