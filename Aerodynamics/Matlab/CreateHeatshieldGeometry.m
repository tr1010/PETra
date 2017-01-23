

% Create heatshield geometry
NoseLength = % radius of nose sphere
RBody = % Radius of overall body
ShoulderHeight % see diagram

nop=100;

nopNose=10;                                  % Number of mesh points for nose

nopFrustum=nop;                             % Number of mesh points for frustum

nopShoulder=10;                              % Number of mesh points for shoulder

r1=linspace(0,NoseLength,nopNose);                          % Define radial vector for nose

r2=linspace(NoseLength,RBody-ShoulderHeight,nopFrustum);    % Define radial vector for frustum

r3=linspace(RBody-ShoulderHeight,RBody,nopShoulder);        % Define radial vector for shoulder

% Radius points
rho=[r1,r2(2:end),r3(2:end)];                                 % Combine to create single vector for full heatshield radius

% Points around
theta=linspace(0,2*pi,length(rho));                           % Define theta angle vector for full heatshield

[theta,rho]=meshgrid(theta,rho);                                % Create mesh grid in r and theta co-ordinates

 

for i=1:length(rho)                                               

    for j=1:length(rho)

        if rho(i,j)<=NoseLength

            z(i,j)=-sqrt(abs(RNose^2-rho(i,j).^2))+RNose;

         elseif rho(i,j)>NoseLength && rho(i,j)<=(FrustumLength+NoseLength)

            z(i,j)=rho(i,j)/tan(Theta)+ConeVertex;

        elseif rho(i,j)<=RBody

            z(i,j)=-sqrt(abs(RShoulder^2-(rho(i,j)-(RBody-RShoulder)).^2))+(ShoulderLength+FrustumHeight+NoseHeight);

        j=j+1;

        end

    end

    i=i+1;

end

[x,y,z]=pol2cart(theta,r,z);            % Transforms cylindrical co-ords to cartesian co-ords

surf(x,y,z)  

 
