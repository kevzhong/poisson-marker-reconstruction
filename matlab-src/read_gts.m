function [xyzv,xyz,verts_of_face,nhat,Atri] = read_gts(fname)

data = readmatrix(fname,'FileType','text');

Nvert = data(1,1);
Nedge = data(1,2);
Nface = data(1,3);
data = data(2:end,:);


xyzv = data(1:Nvert,:);


vert_of_edge = data(Nvert+1:Nvert+Nedge,:);
vert_of_edge = vert_of_edge(:,1:2); 
edges_of_face = data(Nvert+1+Nedge:end,:); %Edges of face

%%%%%%%%%%%%%%  FINISHED READING GTS FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verts_of_face = zeros(Nface,3);
%%%%%%%% FACE_CON_VERTS     COMPUTE VERTICES OF FACES CONNECTIVITY %%%%%%%%
for i = 1:Nface
    %i=1
    %Vertices involved in a given face: 6 total
    edges = [edges_of_face(i,:)];

    %Edges 1 and 2 of face i
    e1 = edges(1);
    e2 = edges(2); 

    %From set_connectivity(), cf. with gts convention
    if vert_of_edge(e1,2) == vert_of_edge(e2,1)
        verts_of_face(i,1) = vert_of_edge(e1,1);
        verts_of_face(i,2) = vert_of_edge(e1,2);
        verts_of_face(i,3) = vert_of_edge(e2,2);
    elseif vert_of_edge(e1,2) == vert_of_edge(e2,2)
        verts_of_face(i,1) = vert_of_edge(e1,1);
        verts_of_face(i,2) = vert_of_edge(e1,2);
        verts_of_face(i,3) = vert_of_edge(e2,1);
    elseif vert_of_edge(e1,1) == vert_of_edge(e2,1)
        verts_of_face(i,1) = vert_of_edge(e1,2);
        verts_of_face(i,2) = vert_of_edge(e1,1);
        verts_of_face(i,3) = vert_of_edge(e2,2);
    else 
        verts_of_face(i,1) = vert_of_edge(e1,2);
        verts_of_face(i,2) = vert_of_edge(e1,1);
        verts_of_face(i,3) = vert_of_edge(e2,1);
    end

end


Atri = compute_area_tri(verts_of_face,xyzv);
xyz = get_triCentroids(verts_of_face,xyzv);


%%%%%%%%%%%%%%     COMPUTE FACE NORMALS     %%%%%%%%%%%%%%%%%%%
nhat = zeros(Nface,3);
for i = 1:Nface
    p1 = xyzv( verts_of_face(i,1) , : );
    p2 = xyzv( verts_of_face(i,2) , : );
    p3 = xyzv( verts_of_face(i,3) , : );

    e1 = p2 - p1;
    e2 = p3 - p1;
    %e3 = p3 - p2;

    nhatb = cross(e1,e2);

    nhatb = nhatb / norm(nhatb);

    nhat(i,:) = nhatb;
end
