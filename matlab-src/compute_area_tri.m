function Atri = compute_area_tri(F,V)
%Compute triangle areas defined by faces F(Nf,3) and vertices V(Nv,3)
%Nf = size(F,1);

% for i = 1:Nf
%     P1 = [ V( F(i,1) , 1 ) , V( F(i,1) , 2 ) , V( F(i,1) , 3 ) ] %[x,y,z]
%     P2 = [ V( F(i,2) , 1 ) , V( F(i,2) , 2 ) , V( F(i,2) , 3 ) ] %[x,y,z]
%     P3 = [ V( F(i,3) , 1 ) , V( F(i,3) , 2 ) , V( F(i,3) , 3 ) ] %[x,y,z]
% 
%     Atri = 1/2*norm(cross(P2-P1,P3-P1));
% end

P1 = [ V( F(:,1) , 1 ) , V( F(:,1) , 2 ) , V( F(:,1) , 3 ) ] ;%[x,y,z]
P2 = [ V( F(:,2) , 1 ) , V( F(:,2) , 2 ) , V( F(:,2) , 3 ) ] ;%[x,y,z]
P3 = [ V( F(:,3) , 1 ) , V( F(:,3) , 2 ) , V( F(:,3) , 3 ) ] ;%[x,y,z]

Atri = 1/2 * vecnorm(cross(P2-P1,P3-P1),2,2) ;


end