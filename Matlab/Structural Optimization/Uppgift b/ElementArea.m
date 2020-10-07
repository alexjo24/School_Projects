function [ A0 ] = ElementArea( ec )
%Compute element area



x1 = ec(1,1);
x2 = ec(1,2);
x3 = ec(1,3);
x4 = ec(1,4);

y1 = ec(1,5);
y2 = ec(1,6);
y3 = ec(1,7);
y4 = ec(1,8);

%Element Area.
A0 =  (1/2)*abs(x1*y2+x2*y3+x3*y4+x4*y1-x2*y1-x3*y2-x4*y3-x1*y4);    

end

