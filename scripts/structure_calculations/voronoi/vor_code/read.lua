s = {}

function docalculation (p)
    local l = s[0]
    local xoffset = s[1]
    local yoffset = s[2]
    local zoffset = s[3]
    p:addpoint(l, xoffset, yoffset, zoffset)
end
