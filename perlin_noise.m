% perlin noise generators 
% There are two perlin noise generators but the second one doesn't seem to work well. 
% perlin_noise is good but super slow!!!!!!!
function s = perlin_noise(s)
% s: A place holder of the matrix for Perlin noise. 
% Return: a matrix of the same size of the input. 

[n, m] = size(s);
i = 0;
w = sqrt(n*m);

while w > 3
    i = i + 1;
    d = interp2(randn(n, m), i-1, 'spline');
    s = s + i * d(1:n, 1:m);
    w = w - ceil(w/2 - 1);
end
s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
end

%% Not very useful implementation of Perlin noise after this point. 
function noise=perlin(values,x,y,z)
if(numel(values)~=512)
    values=randperm(256)-1;
    values=[values values];
end
x=abs(x);
X=bitand(floor(x),255);
x=x-floor(x);
u=fade(x);
A=values(1+X);
noise=linterp(u,grad1d(values(1+X),x),grad1d(values(1+X+1),x-1));
if(nargin>2)
    y=abs(y);
    Y=bitand(floor(y),255);
    y=y-floor(y);
    v=fade(y);
    A=A+Y;
    B=values(1+X+1)+Y;
    noise=linterp(u,linterp(u,grad2d(values(1+A),x,y),grad2d(values(1+B),x-1,y)),linterp(u,grad2d(values(1+A+1),x,y-1),grad2d(values(1+B+1),x-1,y-1)));
end
if(nargin>3)
    z=abs(z);
    Z=bitand(floor(z),255);
    z=z-floor(z);
    w=fade(z);
    AA=values(1+A)+Z;
    AB=values(1+A+1)+Z;
    BA=values(1+B)+Z;
    BB=values(1+B+1)+Z;
    noise=linterp(  w, ...
        linterp(v, ...
        linterp(u, ...
        grad3d(values(1+AA),x,y,z), ...
        grad3d(values(1+BA),x-1,y,z)), ...
        linterp(u, ...
        grad3d(values(1+AB),x,y-1,z), ...
        grad3d(values(1+BB),x-1,y-1,z))), ...
        linterp(v, ...
        linterp(u, ...
        grad3d(values(1+AA+1),x,y,z-1), ...
        grad3d(values(1+BA+1),x-1,y,z-1)), ...
        linterp(u, ...
        grad3d(values(1+AB+1),x,y-1,z-1), ...
        grad3d(values(1+BB+1),x-1,y-1,z-1))));
end
end

function l=linterp(t,a,b)
l=a+t*(b-a);
end

function t=fade(t)
t=6*t^5-15*t^4+10*t^3;
end

function g=grad1d(hash,x)
if(bitand(hash,1))
    g=-x;
else
    g=x;
end
end

function g=grad2d(hash,x,y)
h=bitand(hash,3);
if(bitand(h,2))
    u=-x;
else
    u=x;
end
if(bitand(h,1))
    v=-y;
else
    v=y;
end
g=u+v;
end

function g=grad3d(hash,x,y,z)
h=bitand(hash,15);
if(h<8)
    u=x;
else
    u=y;
end
if(h<4)
    v=y;
elseif(h==12 || h==14)
    v=x;
else
    v=z;
end
if(bitand(h,1))
    if(bitand(h,2))
        g=-u-v;
    else
        g=-u+v;
    end
else
    if(bitand(h,2))
        g=u-v;
    else
        g=u+v;
    end
end
end