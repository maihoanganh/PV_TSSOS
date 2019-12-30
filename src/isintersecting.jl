function isintersecting(v1, v2)
    flag=0
    for v in v1
        if v in v2 
            flag=1
            break
        end
    end
return flag 
end