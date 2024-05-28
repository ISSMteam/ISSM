function rift=updateriftindexing(rift,elconv,nodeconv)
%UPDATERIFTINDEXING - update rift indexing, using mesh to new mesh conversion tables
%     See also meshaddrift

rift.segments(:,1:2)=nodeconv(rift.segments(:,1:2));
rift.segments(:,3)=elconv(rift.segments(:,3));
rift.pairs=elconv(rift.pairs);
rift.tips=nodeconv(rift.tips);

rift.penaltypairs(:,1:2)=nodeconv(rift.penaltypairs(:,1:2));
rift.penaltypairs(:,3:4)=elconv(rift.penaltypairs(:,3:4));
