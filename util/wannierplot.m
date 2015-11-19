function wannierplot(Path,Seedname,Isovalue,Wannieridx)
%plots the wannier orbitals

%initialization
Driver = strcat(Seedname,'_mdriver');
% Path = feval(Driver,'path');
addpath(Path)
% Driver = strcat(Path,Driver);
num_wann = feval(Driver,'num_wann');
ConvBasis = feval(Driver,'basis');%Angstrom
Atom_pos = feval(Driver,'atom_positions');
for i=1:size(Atom_pos,1)
   Atom_pos(i,:) = Atom_pos(i,1)*ConvBasis(:,1)' +...
                   Atom_pos(i,2)*ConvBasis(:,2)' + ...
                   Atom_pos(i,3)*ConvBasis(:,3)' ;
end
Atom_rad = feval(Driver,'atom_radii');
Multiplicities = feval(Driver,'multiplicities');
tmp = [];
Atom_col = [];
coloridx =0;

for i=1:length(Multiplicities)
    tmp = [tmp; repmat(Atom_rad(i),Multiplicities(i),1)];
    coloridx = coloridx + 1;
    Atom_col = [Atom_col;  repmat(get_color(coloridx),Multiplicities(i),1)];
end
Atom_rad = tmp;
PlotOrigin = feval(Driver,'grid_origin');
% Origin = %floor(PlotOrigin)* ConvBasis;
PlotOrigin = PlotOrigin * ConvBasis;
Origin = PlotOrigin;
Gridvec = feval(Driver,'grid_vectors');
plotmode = feval(Driver,'plotmode'); 
printmode = feval(Driver,'printmode');

Filename_psink1 = strcat(Path,Seedname,'_',int2str(Wannieridx),'.psink');
if strcmp(plotmode,'phase')
   Filename_psink1 = strcat(Path,Seedname,'_',int2str(Wannieridx),'.psink');
   Filename_psiarg1 = strcat(Path,Seedname,'_',int2str(Wannieridx),'.psiarg');
end

%read-in number of grid points
fid1=fopen(Filename_psink1);
fgets(fid1);
tmp = str2num(fgets(fid1));
Npoints(1) = tmp(1);
tmp = str2num(fgets(fid1));
Npoints(2) = tmp(1);
tmp = str2num(fgets(fid1));
Npoints(3) = tmp(1);
fclose(fid1);
wien_switch =0;
%read-in data
if wien_switch
    [a,b,c,d,e]=textread(Filename_psink1,'%f%f%f%f%f','headerlines',5);   
    Data.Abs = reshape([a b c d e]',1,5*length(a));
else
    [a,b,c,d,e,f,g,h,i,j]=textread(Filename_psink1,'%f%f%f%f%f%f%f%f%f%f','headerlines',5);   
    Data.Abs = reshape([a b c d e f g h i j]',1,10*length(a));
end
    
if strcmp(plotmode,'phase')
   [a,b,c,d,e,f,g,h,i,j]=textread(Filename_psiarg1,'%f%f%f%f%f%f%f%f%f%f','headerlines',0);
   Data.Arg = reshape([a b c d e f g h i j]',1,10*length(a));
end


%prepare plot
clf

x1 = ConvBasis*Gridvec(:,1);
x2 = ConvBasis*Gridvec(:,2);
x3 = ConvBasis*Gridvec(:,3);

x = PlotOrigin(1):abs(x1(1)-PlotOrigin(1))/(Npoints(1)-1):abs(x1(1));
y = PlotOrigin(2):abs(x2(2)-PlotOrigin(2))/(Npoints(2)-1):abs(x2(2));
z = PlotOrigin(3):abs(x3(3)-PlotOrigin(3))/(Npoints(3)-1):abs(x3(3));

[X,Y,Z] = ndgrid(x,y,z);

Psi.Abs = zeros(Npoints(1),Npoints(2),Npoints(3));
Psi.Arg = zeros(Npoints(1),Npoints(2),Npoints(3));
runner=0;
for i1=1:Npoints(1)
    for i2=1:Npoints(2)
        for i3=1:Npoints(3)
            runner = runner + 1;

                 Psi.Abs(i1,i2,i3) = Data.Abs(runner);
                 if strcmp(plotmode,'phase')
                    Psi.Arg(i1,i2,i3) = Data.Arg(runner);
                 end

              
            
        end
    end
end

nx=floor(max(abs(Gridvec(1,:))));
ny=floor(max(abs(Gridvec(2,:))));
nz=floor(max(abs(Gridvec(3,:))));
direction = diag(sign([x1';x2';x3']-repmat(Origin,3,1)));

for ix=0:nx
    for iy=0:ny
        for iz=0:nz
            
            vec = Origin' + ix*direction(1)*ConvBasis(:,1)+iy*direction(2)*ConvBasis(:,2)+iz*direction(3)*ConvBasis(:,3);
            plot_unitcell(vec,Atom_pos,Atom_rad,Atom_col);
        end
    end
end


col = lines;
col = col(42,:);

if strcmp(plotmode,'phase')
 
   Psi.Arg = cos(Psi.Arg);
   [faces,verts,colors] = isosurface(X,Y,Z,Psi.Abs,Isovalue,Psi.Arg);
   pa=patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,'FaceColor','interp','EdgeColor','none','Facealpha',0.8); 
else
   [faces,verts] = isosurface(X,Y,Z,Psi.Abs,Isovalue);
   pa=patch('Vertices',verts,'Faces',faces,'FaceColor',col,'EdgeColor','none','Facealpha',0.8);
end
camlight; 
lighting phong
colormap(jet)
if strcmp(plotmode,'phase')
  cb = colorbar;
  set(get(cb,'title'),'string','cos(psi)');
  caxis([-1 1])%([min(min(min(Psi.Arg))) max(max(max(Psi.Arg)))])
end


%axis([Origin(1) x1(1) Origin(2) x2(2) Origin(3) x3(3)])
orient landscape
try
   axis1 = feval(Driver,'axis');
   axis2 = [axis1([1 3 5])*ConvBasis axis1([2 4 6])*ConvBasis];
   axis(axis2([1 4 2 5 3 6]))
   
catch
   axis([x(1) x(end) y(1) y(end) z(1) z(end)]);
   axis equal
end

xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
hold off

if ~strcmp(printmode,'off')
  
  Filename_print = strcat(Path,Seedname,'_',int2str(Wannieridx));
  view(3)
  print(strcat(Filename_print,'_plot1.eps'),'-depsc2')
  print(strcat(Filename_print,'_plot1.jpg'),'-djpeg100')
  view(2)
  print(strcat(Filename_print,'_plot2.eps'),'-depsc2')
  print(strcat(Filename_print,'_plot2.jpg'),'-djpeg100')
  view(0,0)
  print(strcat(Filename_print,'_plot3.eps'),'-depsc2')
  print(strcat(Filename_print,'_plot3.jpg'),'-djpeg100')
  view(147,48)
  print(strcat(Filename_print,'_plot4.eps'),'-depsc2')
  print(strcat(Filename_print,'_plot4.jpg'),'-djpeg100')
  view(-10,62)
  print(strcat(Filename_print,'_plot5.eps'),'-depsc2')
  print(strcat(Filename_print,'_plot5.jpg'),'-djpeg100')
end



function plot_unitcell(vec,Atom_pos,Atom_rad,Atom_col)
for i=1:size(Atom_pos,1)
    drawSphere(vec'+Atom_pos(i,:),Atom_rad(i),Atom_col(i,:),'LineStyle','none');
    hold on
end

function col=get_color(idx)
map =lines;
col = map(idx,:);

%--------------------------------------------------------------------------
function varargout = drawSphere(varargin)
%DRAWSPHERE draw a sphere as a mesh
%
%   drawSphere(XC, YC, ZC, R)
%   drawSphere([XC YC ZC], R)
%   drawSphere([XC YC ZC R])
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005
%

%   HISTORY
%   2006-05-19 : use centered sphere with radius 1 when no input specified
%   04/01/2007: typo

options = {};
for i=1:length(varargin)
    if ischar(varargin{i})        
        options = varargin(i:end);
        varargin = varargin(1:i-1);
        break;
    end
end


if isempty(varargin)
    xc = 0;	yc = 0; zc = 0;
    r = 1;    
elseif nargin==1
    sphere = varargin{1};
    xc = sphere(:,1);
    yc = sphere(:,2);
    zc = sphere(:,3);
    r  = sphere(:,4);
    col = 'g';
elseif nargin==2
    center = varargin{1};
    xc = center(1);
    yc = center(2);
    zc = center(3);
    r  = varargin{2};
    col = 'g';
elseif nargin>2
    center = varargin{1};
    xc = center(1);
    yc = center(2);
    zc = center(3);
    r  = varargin{2};
    col= varargin{3};
   
else
    error('drawSphere : please specify center and radius');
end

nphi = 100;
ntheta = 100;

theta = (0:ntheta)/ntheta*pi;
phi = (0:nphi)/nphi*2*pi;

sintheta = sin(theta);

x = xc + cos(phi')*sintheta*r;
y = yc + sin(phi')*sintheta*r;
z = zc + ones(length(phi),1)*cos(theta)*r;


if nargout == 0
    surf(x,y,z, 'FaceColor',col, options{:}); 
elseif nargout == 1
    varargout{1} = surf(x,y,z, 'FaceColor', col, options{:});
elseif nargout == 3
    varargout{1} = x; 
    varargout{2} = y; 
    varargout{3} = z; 
end




