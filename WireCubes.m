function [ model_handle ] = WireCubes( center,cubesize,color,Alpha )
%WireCubes plots cubic wireframe according to the cube cener and its size
%INPUT: 
%size: nX1 vector containing the size of each cubic element
%centers: nX3 matrix containing the x,y,z coordinates of the cubes
%OUTPUT:
%model_handle: an object handle to the plotted model
NumberOfCubes=size(center,1);
verticesTemplate=[0 0 0;
                                0 1 0;
                                1 1 0;
                                1 0 0;
                                0 0 1;
                                0 1 1;
                                1 1 1;
                                1 0 1];
facesTemplate=[1 2 3 4;
                        5 6 7 8;
                         3 4 8 7;
                         1 2 6 5;
                         2 3 7 6;
                          1 4 8 5];
 FV.faces=zeros(NumberOfCubes*6,4);  
 FV.vertices=zeros(NumberOfCubes*8,3);  
 faceindex=1;
 vertexindex=1;
for i=1:NumberOfCubes
    verts = (verticesTemplate-0.5).*cubesize(i,1)+repmat(center(i,:),8,1);
    faces=facesTemplate+(i-1)*8;
      FV.vertices(vertexindex:vertexindex+7,1:3)=verts;
      FV.faces(faceindex:faceindex+5,1:4)=faces;
      vertexindex=vertexindex+8;
      faceindex=faceindex+6;
  end
model_handle=patch(FV,'edgecolor',color,'facecolor',color);
alpha(Alpha);
end
