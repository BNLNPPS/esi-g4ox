import os

# Define the file path
file_path = "/esi/eic-opticks/CSG/csg_intersect_leaf_zsphere.h"

# Define the old and new code blocks
old_code = '''\
if(valid_isect)
    {
        isect.w = t_cand ;
        if( t_cand == t1sph || t_cand == t2sph)
        {
            isect.x = (O.x + t_cand*D.x)/radius ; // normalized by construction
            isect.y = (O.y + t_cand*D.y)/radius ;
            isect.z = (O.z + t_cand*D.z)/radius ;
        }
        else
        {
            isect.x = 0.f ;
            isect.y = 0.f ;
            isect.z = t_cand == t_PCAP ? -1.f : 1.f ;
        }
    }'''

new_code = '''\
if(valid_isect)
{
    isect.w = t_cand ;
    if( t_cand == t1sph || t_cand == t2sph)
    {
        // Reflect from inner spherical surface: flip normal
        isect.x = -(O.x + t_cand*D.x)/radius ;
        isect.y = -(O.y + t_cand*D.y)/radius ;
        isect.z = -(O.z + t_cand*D.z)/radius ;
    }
    else
    {
        // Reflect from inner caps: flip Z normal
        isect.x = 0.f ;
        isect.y = 0.f ;
        isect.z = t_cand == t_PCAP ? 1.f : -1.f ;
    }
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
        print("ZSphere intersection code replaced successfully.")
    else:
        print("DID NOT FIND CONE INTERSECTION TO BE REPLACED")
