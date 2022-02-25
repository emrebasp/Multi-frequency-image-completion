function result = normalizer(array)

if length(size(array)) == 2
mx = max(max(array)); 
mn = min(min(array));
end

if length(size(array)) == 3
mx = max(max(max(array))); 
mn = min(min(min(array)));
end

if length(size(array)) == 4
mx = max(max(max(max(array)))); 
mn = min(min(min(min(array))));
end

if length(size(array)) == 5
mx = max(max(max(max(max(array))))); 
mn = min(min(min(min(min(array)))));
end



result = (array - mn)/(mx-mn);

end