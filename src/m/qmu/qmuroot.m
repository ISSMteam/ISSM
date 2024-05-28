function root=qmuroot(string)
%QMUROOT - return root of a distributed descriptor

root='';
found=0;
for i=1:length(string),
	if ((49<=double(string(i))) && (double(string(i)<=57)))
		break;
	else
		root=[root string(i)];
	end
end
