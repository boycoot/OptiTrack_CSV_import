function data=optitrack_import(name,body_names)

% Loads data from OptiTrack CSV files
%
% name: CSV filename (string)
% body_names: body names that should be extracted from the log (cell of strings)
%
% Example:
% name='Take 2016-07-14 08.56.30 PM V2p0 PIV second shot';
% body_names={'DelFly_rigid','DelFly_flexible','Synchronisation'};
%
% The CSV file should be exported with default settings, e.g. y points up and quaternions are used for orientation.
% The output has z axis pointing up and uses standard airspace roll-pitch-yaw convention
%
% 
% Author: Matej Karasek (matejkarasek@gmail.com), 2016

Nbodies=length(body_names);

openedFile = fopen([name '.csv'],'r');

% read the first 7 lines and save individiual cells delimited by comma
for i=1:7
    line=fgetl(openedFile);
    lines(i).string=regexp(line, ',', 'split');
end

% parse the general properties
data.formatVersion=str2double(lines(1).string(2));
data.fps=str2double(lines(1).string(8));
data.date=lines(1).string(10);
data.Nframes=str2double(lines(1).string(14));
data.filename=name;

if data.formatVersion~=1.21
    disp('Warning: this format version of OptiTrack csv file has not been tested, please check the results carefully')
end

% create index of all the columns
index=1:length(lines(4).string);

% load the rest of the data
data.CSVdata=csvread([name '.csv'],7,0);
data.frames=data.CSVdata(:,1);
data.time=data.CSVdata(:,2);

% for each body, find the relevant data
for n=1:Nbodies
    data.body(n).name=body_names(n);
    Nchar=length(char(body_names(n)));
    
    % for each body, find the relevant columns
    nameMatchFull=strcmp(lines(4).string,char(body_names(n)));
    nameMatchPartial=strncmp(lines(4).string,char(body_names(n)),Nchar);
    markerMatch=strcmp(lines(3).string,'Marker');
    rigidBodyMarkerMatch=strcmp(lines(3).string,'Rigid Body Marker');
    
    % body columns
    data.body(n).cBody=index(nameMatchFull);
    % rigid body (fitted) marker columns (4 columns per marker)
    data.body(n).cMarkersFitted=index((nameMatchPartial-nameMatchFull).*rigidBodyMarkerMatch==1);
    Nmarkers=length(data.body(n).cMarkersFitted)/4; % 4 columns per marker
    % raw marker columns (3 columns per marker)
    data.body(n).cMarkers=index((nameMatchPartial-nameMatchFull).*markerMatch==1);
    
    % set the lines when body was untracked (marker error==0) to NaN
    nanframes=data.frames(data.CSVdata(:,data.body(n).cBody(8))==0)+1;
    data.CSVdata(nanframes,index(nameMatchPartial==1))=NaN;
    
    % body position, orientation and mean marker error
    data.body(n).pos=data.CSVdata(:,data.body(n).cBody([7 5 6])); % OptiTrack z,x,y --> x,y,z
    data.body(n).quat=data.CSVdata(:,data.body(n).cBody(1:4));
    data.body(n).meanError=data.CSVdata(:,data.body(n).cBody(8));
    
    % marker positions
    j=0;
    for i=1:Nmarkers
        data.body(n).marker(i).posFitted=data.CSVdata(:,data.body(n).cMarkersFitted([4*i-1 4*i-3 4*i-2])); % OptiTrack z,x,y --> x,y,z
        data.body(n).marker(i).qual=data.CSVdata(:,data.body(n).cMarkersFitted(4*i));
        
        if sum(abs(data.body(n).marker(i).qual))==0 % the matching marker was not seen in the entire recording
            data.body(n).marker(i).pos=NaN(data.Nframes,3);
        else
            j=j+1; % increment only if the marker was seen
            data.body(n).marker(i).pos=data.CSVdata(:,data.body(n).cMarkers([3*j 3*j-2 3*j-1])); % OptiTrack z,x,y --> x,y,z
        end
    end
    
    % calculate roll, pitch and yaw
    roll=NaN(size(data.time));
    pitch=NaN(size(data.time));
    yaw=NaN(size(data.time));
    
    iprev=1;
    
    qx=data.body(n).quat(:,1);
    qy=data.body(n).quat(:,2);
    qz=data.body(n).quat(:,3);
    qw=data.body(n).quat(:,4);
    
    for i=1:data.Nframes
        % rotation matrix from OptiTrack quaternion
        R=[1-2*(qy(i)^2+qz(i)^2)        2*(qx(i)*qy(i)-qz(i)*qw(i))  2*(qx(i)*qz(i)+qy(i)*qw(i))
            2*(qx(i)*qy(i)+qz(i)*qw(i))  1-2*(qx(i)^2+qz(i)^2)        2*(qy(i)*qz(i)-qx(i)*qw(i))
            2*(qx(i)*qz(i)-qy(i)*qw(i))  2*(qy(i)*qz(i)+qx(i)*qw(i))  1-2*(qx(i)^2+qy(i)^2)      ];
        
        % orientation from OptiTrack log (roll around OptiTrack Z, pitch around OptiTrack X, Yaw around OptiTrack Y)
        roll(i,1)=atan2d(R(2,1),R(2,2));
        pitch(i,1)=atan2d(-R(2,3),real(sqrt(1-R(2,3)^2))); % real added to avoid complex numbers (most likely due to rounding errors)
        yaw(i,1)=atan2d(R(1,3),R(3,3));
        
        % making yaw continuous
        if ~isnan(yaw(i,1))
            if i>1
                while abs(yaw(i,1)-yaw(iprev,1))>180
                    if yaw(i,1)-yaw(iprev,1)>180
                        yaw(i,1)=yaw(i,1)-360;
                    elseif yaw(i,1)-yaw(iprev,1)<-180
                        yaw(i,1)=yaw(i,1)+360;
                    end
                end
                iprev=i;
            end
        end
        
        RE2B=Rot_x(-roll(i)/180*pi)*Rot_y(-pitch(i)/180*pi)*Rot_z(-yaw(i)/180*pi);
        
        % transform marker positions to body axes
        for jj=1:Nmarkers
            data.body(n).marker(jj).posBody(i,:)=(RE2B*(data.body(n).marker(jj).pos(i,:)-data.body(n).pos(i,:))')';
            data.body(n).marker(jj).posBodyFitted(i,:)=(RE2B*(data.body(n).marker(jj).posFitted(i,:)-data.body(n).pos(i,:))')';
        end
        
        if i/1000==round(i/1000)
            disp(['body ' num2str(n) ', frame ' num2str(i)])
        end
    end
    
    % store the roll, pitch and yaw angles in the structure
    data.body(n).roll=roll;
    data.body(n).pitch=pitch;
    data.body(n).yaw=yaw;
end
