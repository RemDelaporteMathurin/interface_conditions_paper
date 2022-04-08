%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TMAP RESULTS CONVERSION
%%%%%%%% 2020.04.23 LSPM JM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
warning off
format long

% RUN WITH MATLAB or OCTAVE (free tools WIN/MAC/LINUX)
% https://www.gnu.org/software/octave/download.html

%%%%%%%% MUST BE COMPLETED %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of the file
ftoread="./extract.txt";

% Nomber of segments
nb_segment=3;

% Please insert the number of cell and size for each segments

s.nb_cell_1=[8,8,4,8,4,8,4,8,4,29];
s.size_cell_1 = [2.5e-9, 1e-8,   2.5e-8, 1e-7, 2.5e-7,...
                1e-6, 2.5e-6, 1e-5, 2.5e-5, 2e-4];
          
s.nb_cell_2=[100];
s.size_cell_2 = [1e-5];

s.nb_cell_3=[100];
s.size_cell_3 = [1.5e-5];

% EXAMPLE :
% FROM TMAP 
%        delx = 0.0,8*2.5e-9,8*1e-8,4*2.5e-8,8*1e-7,4*2.5e-7,8*1e-6                        
%              4*2.5e-6,8*1e-5,4*2.5e-5,29*2e-4,0.0,end
% PUT HERE 
%         s.nb_cell_i = [8,8,4,8,4,8,4,8,4,29];
%         s.size_cell_i = [2.5e-9, 1e-8,   2.5e-8, 1e-7,   2.5e-7,...
%                            1e-6, 1e-5,   2.5e-6, 2.5e-5, 2e-4];
% where i is the number of the segment
                            
% Number of trap(s) per segments
traps_per_segment=[2,1,2];


%%%%%%%% NOT NEED TO BE EDITED BELOW %%%%%%%%%%%%%%%%%

nb_cells_per_s=[sum(s.nb_cell_1)+2];
for i=2:nb_segment
  nb_cells_per_s=[ nb_cells_per_s sum(s.(['nb_cell_',num2str(i)]))+2];
end

disp('Number of nodes per segment :');
nb_cells_per_s
disp(['Total number of nodes = ',num2str(sum(nb_cells_per_s))]);

dx=ones(sum(nb_cells_per_s),1);
k=1;
for m=1:nb_segment
dx(k)=0;
k=k+1;
current=s.(['nb_cell_',num2str(m)]);
current2=s.(['size_cell_',num2str(m)]);
for i=1:length(current)
    for j=1:current(i)
        dx(k)=current2(i);
        k=k+1;
    end
end
dx(k)=0;
k=k+1;
end
dx(1:90);

x=zeros(sum(nb_cells_per_s),1);
x(1)=0;
for i=2:sum(nb_cells_per_s)
  x(i)=x(i-1)+(dx(i-1)+dx(i))/2;
end

disp(['Segment 1: [',num2str(x(1)),':',num2str(x(nb_cells_per_s(1))),']']);
for i=2:nb_segment
disp(['Segment ',num2str(i)',': [',num2str(x(sum(nb_cells_per_s(1:i-1)))),...
                               ':',num2str(x(sum(nb_cells_per_s(1:i)))),']']);
end



temperature=zeros(sum(nb_cells_per_s),1);
mobile=zeros(sum(nb_cells_per_s),1);
trapped=zeros(sum(nb_cells_per_s),sum(traps_per_segment));
mobile_inventory=zeros(numel(nb_cells_per_s),1);
trapped_inventory=zeros(numel(nb_cells_per_s),1);

fid = fopen(ftoread);
time=str2num(fgetl(fid)(15:25));

m=1;
for k=1:nb_segment
fgetl(fid);
fgetl(fid);
fgetl(fid);
for i=1:nb_cells_per_s(k)
if(k==1)
    temperature(i)=fscanf(fid, "%f",1);
else
    temperature(i+sum(nb_cells_per_s(1:k-1)))=fscanf(fid, "%f",1);
end
end
fgetl(fid);
fgetl(fid);
fgetl(fid);
for i=1:nb_cells_per_s(k)
if(k==1)
    mobile(i)=fscanf(fid, "%f",1);
else
    mobile(i+sum(nb_cells_per_s(1:k-1)))=fscanf(fid, "%f",1);  
end    
end
fgetl(fid);
mobile_inventory(k)=str2num(fgetl(fid)(26:36));
fgetl(fid);
fgetl(fid);
for t=1:traps_per_segment(k)
  fgetl(fid);
  for i=1:nb_cells_per_s(k)
     if(k==1)
        trapped(i,m)=fscanf(fid, "%f",1);
     else
        trapped(i+sum(nb_cells_per_s(1:k-1)),m)=fscanf(fid, "%f",1);
     end
  end
  m++;
  fgetl(fid);
end  
trapped_inventory(k)=str2num(fgetl(fid)(53:63));
end
fclose (fid);

disp('Inventory per segment (at/m2) :');
mobile_inventory
trapped_inventory


fid = fopen(['convert-',num2str(time),'.txt'],'w');
fprintf(fid,"node\tx (m)    \tTemp. (K) \tMobile (at/m3)")
  for j=1:sum(traps_per_segment)
     fprintf(fid,"\tTrap %d (at/m3)  ",j);
  end
fprintf(fid,"\n");
for i=1:length(x)
  fprintf(fid,"%d\t%E\t%E\t%E",i,x(i),temperature(i),mobile(i));
  for j=1:sum(traps_per_segment)
     fprintf(fid,"\t%E",trapped(i,j));
  end
  fprintf(fid,"\n");
end
fclose (fid);

figure(1)
hold on
grid on
plot(x,temperature)
xlabel('x (m)')
ylabel('Temperature (K)')
print(['Temperature-',num2str(time),'s.png'], '-dpng')


figure(2)
hold on
grid on
semilogy(x,mobile)
xlabel('x (m)')
ylabel('Mobile concentration (at/m3)')
print(['Mobile-',num2str(time),'s.png'], '-dpng')

figure(3)
hold on
grid on
semilogy(x,trapped)
xlabel('x (m)')
ylabel('Trapped concentration (at/m3)')
print(['Trapped-',num2str(time),'s.png'], '-dpng')




