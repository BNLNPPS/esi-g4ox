import os

# Define the file path
file_path = "/esi/eic-opticks/CSG/csg_intersect_leaf_newcone.h"

# Define the old and new code blocks
old_code = '''\
    valid_isect = t_cand > t_min && t_cand < RT_DEFAULT_MAX ;
    if(valid_isect)
    {
        if( t_cand == t_cap1 || t_cand == t_cap2 )
        {
            isect.x = 0.f ; 
            isect.y = 0.f ;
            isect.z = t_cand == t_cap2 ? 1.f : -1.f  ;   
        }
        else
        { 
            float3 n = normalize(make_float3( o.x+t_cand*d.x, o.y+t_cand*d.y, (z0-(o.z+t_cand*d.z))*tth2  ))  ; 
            isect.x = n.x ; 
            isect.y = n.y ;
            isect.z = n.z ; 
        }
        isect.w = t_cand ; 
    }'''

new_code = '''\
    valid_isect = (t_cand > t_min && t_cand < RT_DEFAULT_MAX);
    if (valid_isect)
    {
        float3 intersection_point = make_float3(
            o.x + t_cand * d.x,
            o.y + t_cand * d.y,
            o.z + t_cand * d.z
        );

        printf("// Intersection at t = %.4f, position: (%.4f, %.4f, %.4f)\\n",
            t_cand,
            intersection_point.x,
            intersection_point.y,
            intersection_point.z
        );

        float3 n = normalize(make_float3(
            intersection_point.x,
            intersection_point.y,
            (z0 - intersection_point.z)*tth2
        ));

        printf("// Intersection with cone side. Normal: (%.4f, %.4f, %.4f)\\n",
            n.x, n.y, n.z
        );

        isect.x = n.x;
        isect.y = n.y;
        isect.z = n.z;
        isect.w = t_cand;
    }'''

# Read the file and replace the target code
if not os.path.exists(file_path):
    print("ERROR: File does not exist:", file_path)
else:
    with open(file_path, 'r') as file:
        content = file.read()

    if old_code in content:
        content = content.replace(old_code, new_code)
        with open(file_path, 'w') as file:
            file.write(content)
        print("Cone intersection code replaced successfully.")
    else:
        print("DID NOT FIND CONE INTERSECTION TO BE REPLACED")
