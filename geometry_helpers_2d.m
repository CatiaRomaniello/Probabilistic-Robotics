1;
function v=t2v(A)
	v(1:2, 1)=A(1:2,3);
	v(3,1)=atan2(A(2,1),A(1,1));
end

% computes the homogeneous transform matrix A of the pose vector v
% A:[ R t ] 3x3 homogeneous transformation matrix, r translation vector
% v: [x,y,theta]  2D pose vector
function A=v2t(v)
  	c=cos(v(3));
  	s=sin(v(3));
	A=[c, -s, v(1) ;
	s,  c, v(2) ;
	0   0  1  ];
end

function A=v2t_3d(v)
  	c=cos(v(6));
  	s=sin(v(6));
	A=[c, -s, 0, v(1) ;
	s,  c, 0, v(2) ;
	0, 0, 1, 0;
	0,   0, 0, 1 ];
end

function v=flattenIsometryByColumns_2d(T)
	v=zeros(6,1);
	v(1:4)=reshape(T(1:2,1:2),4,1);
	v(5:6)=T(1:2,3);
endfunction
% normalizes and angle between -pi and pi
% th: input angle
% o: output angle
function o = normalizeAngle(th)
	o = atan2(sin(th),cos(th));
end

#computes the derivative of atan2(p.y/p.x);
function J = J_atan2(p)
  n2=transpose(p)*p;
  J= 1./n2 * [-p(2) p(1)];
endfunction
