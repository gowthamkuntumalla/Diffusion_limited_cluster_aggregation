function import_time(numfiles)
%Import time data
%numfiles=418;
time=zeros(numfiles,1);
for i=1:numfiles
    filename=sprintf('data_time%d',i);    
    time(i)=importdata(filename); 
end
save('data_time');
end