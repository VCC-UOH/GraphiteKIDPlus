// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <OGF/CurveKID/algo/mesh_geometry.h>

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <random>
#include <cassert>

namespace OGF {
namespace CurveKIDAlgo {

// =======================================================================
// Internal helpers
// =======================================================================

static GEO::vec3 vertex_point(const GEO::Mesh& mesh, GEO::index_t v) {
    const double* p = mesh.vertices.point_ptr(v);
    return GEO::vec3(p[0], p[1], p[2]);
}

static double triangle_area_3pts(
    const GEO::vec3& p0, const GEO::vec3& p1, const GEO::vec3& p2)
{
    return 0.5 * GEO::length(GEO::cross(p1 - p0, p2 - p0));
}

static double cotan_angle(
    const GEO::vec3& v0, const GEO::vec3& v1, const GEO::vec3& v2)
{
    GEO::vec3 u = v0 - v2;
    GEO::vec3 v = v1 - v2;
    double dot      = GEO::dot(u, v);
    double cross_len = GEO::length(GEO::cross(u, v));
    if (cross_len < 1e-15) return 0.0;
    return dot / cross_len;
}

// =======================================================================
// compute_mesh_stats
// =======================================================================

bool compute_mesh_stats(const GEO::Mesh& mesh, MeshStats& stats) {
    stats = MeshStats{};
    stats.nb_vertices = mesh.vertices.nb();
    stats.nb_facets   = mesh.facets.nb();
    if (mesh.facets.nb() == 0) return false;

    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
        if (mesh.facets.nb_vertices(f) == 3) ++stats.nb_triangles;
        else ++stats.nb_non_tri;
    }
    if (stats.nb_triangles == 0) return false;

    struct PairHash {
        size_t operator()(std::pair<GEO::index_t,GEO::index_t> p) const {
            return std::hash<uint64_t>()(
                (uint64_t)p.first << 32 | (uint64_t)p.second);
        }
    };
    std::unordered_map<std::pair<GEO::index_t,GEO::index_t>,int,PairHash> edge_map;

    double sum_edge=0, sum_area=0;
    double min_e=1e30, max_e=-1e30;
    double min_a=1e30, max_a=-1e30;
    int edge_count=0, border_count=0;

    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
        if (mesh.facets.nb_vertices(f) != 3) continue;
        GEO::index_t v0=mesh.facets.vertex(f,0);
        GEO::index_t v1=mesh.facets.vertex(f,1);
        GEO::index_t v2=mesh.facets.vertex(f,2);
        GEO::vec3 p0=vertex_point(mesh,v0);
        GEO::vec3 p1=vertex_point(mesh,v1);
        GEO::vec3 p2=vertex_point(mesh,v2);

        double area = triangle_area_3pts(p0,p1,p2);
        sum_area += area;
        min_a = std::min(min_a, area);
        max_a = std::max(max_a, area);

        GEO::index_t vs[3]={v0,v1,v2};
        for (int e=0;e<3;++e) {
            GEO::index_t ea=vs[e], eb=vs[(e+1)%3];
            if (ea>eb) std::swap(ea,eb);
            auto key=std::make_pair(ea,eb);
            auto it=edge_map.find(key);
            if (it==edge_map.end()) {
                edge_map[key]=1;
                double len=GEO::length(
                    vertex_point(mesh,eb)-vertex_point(mesh,ea));
                sum_edge+=len;
                min_e=std::min(min_e,len);
                max_e=std::max(max_e,len);
                ++edge_count;
            } else { ++it->second; }
        }
    }
    for (auto& kv:edge_map) if(kv.second==1) ++border_count;

    stats.nb_edges         = (GEO::index_t)edge_count;
    stats.nb_border_edges  = (GEO::index_t)border_count;
    stats.min_edge_length  = min_e;
    stats.max_edge_length  = max_e;
    stats.mean_edge_length = edge_count>0 ? sum_edge/edge_count : 0.0;
    stats.min_face_area    = min_a;
    stats.max_face_area    = max_a;
    stats.mean_face_area   = stats.nb_triangles>0
                             ? sum_area/stats.nb_triangles : 0.0;
    stats.total_surface_area = sum_area;
    return true;
}

void compute_face_normals(
    const GEO::Mesh& mesh, std::vector<GEO::vec3>& normals)
{
    normals.resize(mesh.facets.nb());
    for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
        if (mesh.facets.nb_vertices(f)!=3){normals[f]=GEO::vec3(0,0,1);continue;}
        GEO::vec3 p0=vertex_point(mesh,mesh.facets.vertex(f,0));
        GEO::vec3 p1=vertex_point(mesh,mesh.facets.vertex(f,1));
        GEO::vec3 p2=vertex_point(mesh,mesh.facets.vertex(f,2));
        GEO::vec3 n=GEO::cross(p1-p0,p2-p0);
        double len=GEO::length(n);
        normals[f]=len>1e-15?n/len:GEO::vec3(0,0,1);
    }
}

void compute_face_areas(
    const GEO::Mesh& mesh, std::vector<double>& areas)
{
    areas.resize(mesh.facets.nb());
    for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
        if (mesh.facets.nb_vertices(f)!=3){areas[f]=0;continue;}
        GEO::vec3 p0=vertex_point(mesh,mesh.facets.vertex(f,0));
        GEO::vec3 p1=vertex_point(mesh,mesh.facets.vertex(f,1));
        GEO::vec3 p2=vertex_point(mesh,mesh.facets.vertex(f,2));
        areas[f]=triangle_area_3pts(p0,p1,p2);
    }
}

void enumerate_triangle_edges(
    const GEO::Mesh& mesh, std::vector<EdgeInfo>& edges)
{
    struct PairHash {
        size_t operator()(std::pair<GEO::index_t,GEO::index_t> p) const {
            return std::hash<uint64_t>()(
                (uint64_t)p.first<<32|(uint64_t)p.second);
        }
    };
    std::unordered_map<std::pair<GEO::index_t,GEO::index_t>,int,PairHash> em;
    for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
        if (mesh.facets.nb_vertices(f)!=3) continue;
        for (int e=0;e<3;++e) {
            GEO::index_t a=mesh.facets.vertex(f,e);
            GEO::index_t b=mesh.facets.vertex(f,(e+1)%3);
            if (a>b) std::swap(a,b);
            em[{a,b}]++;
        }
    }
    edges.clear(); edges.reserve(em.size());
    for (auto& kv:em) {
        EdgeInfo ei;
        ei.v0=kv.first.first; ei.v1=kv.first.second;
        ei.length=GEO::length(vertex_point(mesh,ei.v1)-vertex_point(mesh,ei.v0));
        ei.is_boundary=(kv.second==1);
        edges.push_back(ei);
    }
}

// =======================================================================
// STAGE 1 — Curvature estimation
// =======================================================================

void compute_vertex_curvatures(
    const GEO::Mesh& mesh,
    std::vector<VertexCurvature>& vcurv)
{
    GEO::index_t nv = mesh.vertices.nb();
    vcurv.assign(nv, VertexCurvature{});

    std::vector<double> mixed_area(nv, 0.0);
    std::vector<double> H_num(nv, 0.0);
    std::vector<double> angle_sum(nv, 0.0);

    for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
        if (mesh.facets.nb_vertices(f)!=3) continue;
        GEO::index_t vi[3];
        GEO::vec3    p[3];
        for (int k=0;k<3;++k) {
            vi[k]=mesh.facets.vertex(f,k);
            p[k]=vertex_point(mesh,vi[k]);
        }
        double area=triangle_area_3pts(p[0],p[1],p[2]);
        if (area<1e-20) continue;

        double angle[3];
        for (int k=0;k<3;++k) {
            GEO::vec3 u=p[(k+1)%3]-p[k];
            GEO::vec3 v=p[(k+2)%3]-p[k];
            double lu=GEO::length(u),lv=GEO::length(v);
            if (lu<1e-15||lv<1e-15){angle[k]=0.0;continue;}
            double c=GEO::dot(u,v)/(lu*lv);
            c=std::max(-1.0,std::min(1.0,c));
            angle[k]=std::acos(c);
        }
        bool obtuse=(angle[0]>M_PI/2||angle[1]>M_PI/2||angle[2]>M_PI/2);

        for (int k=0;k<3;++k) {
            double ma;
            if (obtuse) {
                ma=(angle[k]>M_PI/2)?area*0.5:area*0.25;
            } else {
                double cot1=cotan_angle(p[k],p[(k+2)%3],p[(k+1)%3]);
                double cot2=cotan_angle(p[k],p[(k+1)%3],p[(k+2)%3]);
                double e1sq=GEO::dot(p[(k+1)%3]-p[k],p[(k+1)%3]-p[k]);
                double e2sq=GEO::dot(p[(k+2)%3]-p[k],p[(k+2)%3]-p[k]);
                ma=(cot1*e2sq+cot2*e1sq)/8.0;
                if (ma<0) ma=area/3.0;
            }
            mixed_area[vi[k]]+=ma;
            angle_sum[vi[k]] +=angle[k];
        }

        for (int k=0;k<3;++k) {
            int i=k,j=(k+1)%3,opp=(k+2)%3;
            double cot_opp=cotan_angle(p[i],p[j],p[opp]);
            GEO::vec3 e=p[j]-p[i];
            double contrib=cot_opp*GEO::dot(e,e);
            H_num[vi[i]]+=contrib;
            H_num[vi[j]]+=contrib;
        }
    }

    for (GEO::index_t v=0;v<nv;++v) {
        double A=mixed_area[v];
        if (A<1e-20) continue;
        double K=(2.0*M_PI-angle_sum[v])/A;
        double H=std::abs(H_num[v]/(4.0*A));
        double disc=H*H-K;
        double kappa1,kappa2;
        if (disc<0.0) { kappa1=kappa2=H; }
        else { double sq=std::sqrt(disc); kappa1=H-sq; kappa2=H+sq; }
        vcurv[v].gaussian=K;
        vcurv[v].mean=H;
        vcurv[v].kappa1=kappa1;
        vcurv[v].kappa2=kappa2;
    }
}

void compute_face_curvatures(
    const GEO::Mesh& mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<FaceCurvature>& fcurv)
{
    fcurv.resize(mesh.facets.nb());
    for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
        if (mesh.facets.nb_vertices(f)!=3) continue;
        FaceCurvature fc{};
        for (int k=0;k<3;++k) {
            GEO::index_t v=mesh.facets.vertex(f,k);
            fc.gaussian+=vcurv[v].gaussian;
            fc.mean    +=vcurv[v].mean;
            fc.kappa1  +=vcurv[v].kappa1;
            fc.kappa2  +=vcurv[v].kappa2;
        }
        fc.gaussian/=3.0; fc.mean/=3.0;
        fc.kappa1/=3.0;   fc.kappa2/=3.0;
        fcurv[f]=fc;
    }
}

// Least-squares sphere fit
static bool fit_sphere_ls(
    const std::vector<GEO::vec3>& pts,
    GEO::vec3& center, double& radius, double& rms)
{
    int n=(int)pts.size();
    if (n<4) return false;
    double M[4][4]={};
    double rhs[4]={};
    for (auto& p:pts) {
        double row[4]={2*p.x,2*p.y,2*p.z,-1.0};
        double bi=p.x*p.x+p.y*p.y+p.z*p.z;
        for (int i=0;i<4;++i) {
            rhs[i]+=row[i]*bi;
            for (int j=0;j<4;++j) M[i][j]+=row[i]*row[j];
        }
    }
    double A[4][5];
    for (int i=0;i<4;++i){for(int j=0;j<4;++j)A[i][j]=M[i][j];A[i][4]=rhs[i];}
    for (int col=0;col<4;++col) {
        int pivot=col;
        for (int row=col+1;row<4;++row)
            if(std::abs(A[row][col])>std::abs(A[pivot][col])) pivot=row;
        if(std::abs(A[pivot][col])<1e-15) return false;
        std::swap(A[col],A[pivot]);
        double inv=1.0/A[col][col];
        for (int row=0;row<4;++row) {
            if(row==col) continue;
            double factor=A[row][col]*inv;
            for(int k=col;k<=4;++k) A[row][k]-=factor*A[col][k];
        }
        for(int k=col+1;k<=4;++k) A[col][k]*=inv;
    }
    double cx=A[0][4],cy=A[1][4],cz=A[2][4],d=A[3][4];
    center=GEO::vec3(cx,cy,cz);
    double r2=cx*cx+cy*cy+cz*cz-d;
    if(r2<0) return false;
    radius=std::sqrt(r2);
    double sum=0.0;
    for(auto& p:pts){double err=GEO::length(p-center)-radius;sum+=err*err;}
    rms=std::sqrt(sum/n);
    return true;
}

bool fit_sphere_to_neighborhood(
    const GEO::Mesh& mesh, GEO::index_t f,
    GEO::vec3& center, double& radius, double& residual)
{
    std::unordered_set<GEO::index_t> vset;
    for(int k=0;k<3;++k) vset.insert(mesh.facets.vertex(f,k));
    for(GEO::index_t ff=0;ff<mesh.facets.nb();++ff) {
        if(mesh.facets.nb_vertices(ff)!=3) continue;
        for(int k=0;k<3;++k) {
            if(vset.count(mesh.facets.vertex(ff,k))) {
                for(int j=0;j<3;++j)
                    vset.insert(mesh.facets.vertex(ff,j));
                break;
            }
        }
    }
    std::vector<GEO::vec3> pts;
    pts.reserve(vset.size());
    for(GEO::index_t v:vset) pts.push_back(vertex_point(mesh,v));
    return fit_sphere_ls(pts,center,radius,residual);
}

// =======================================================================
// STAGE 2 — Patch growing
// =======================================================================

GEO::index_t grow_patches(
    const GEO::Mesh& mesh,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t target_count,
    std::vector<GEO::index_t>& patch_id)
{
    GEO::index_t nf=mesh.facets.nb();
    patch_id.assign(nf,GEO::index_t(-1));

    // Build face adjacency
    struct PairHash {
        size_t operator()(std::pair<GEO::index_t,GEO::index_t> p) const {
            return std::hash<uint64_t>()(
                (uint64_t)p.first<<32|(uint64_t)p.second);
        }
    };
    std::unordered_map<
        std::pair<GEO::index_t,GEO::index_t>,
        std::vector<GEO::index_t>,PairHash> edge_to_faces;

    for(GEO::index_t f=0;f<nf;++f) {
        if(mesh.facets.nb_vertices(f)!=3) continue;
        for(int e=0;e<3;++e) {
            GEO::index_t a=mesh.facets.vertex(f,e);
            GEO::index_t b=mesh.facets.vertex(f,(e+1)%3);
            if(a>b) std::swap(a,b);
            edge_to_faces[{a,b}].push_back(f);
        }
    }

    std::vector<std::vector<GEO::index_t>> fadj(nf);
    for(auto& kv:edge_to_faces) {
        auto& fl=kv.second;
        for(GEO::index_t i=0;i<fl.size();++i)
            for(GEO::index_t j=i+1;j<fl.size();++j) {
                fadj[fl[i]].push_back(fl[j]);
                fadj[fl[j]].push_back(fl[i]);
            }
    }

    // Count triangles
    GEO::index_t nb_tri=0;
    for(GEO::index_t f=0;f<nf;++f)
        if(mesh.facets.nb_vertices(f)==3) ++nb_tri;
    GEO::index_t actual_count=std::min(target_count,nb_tri);

    // Furthest-point seeding in curvature space
    double k1_max=1e-10,k2_max=1e-10;
    for(GEO::index_t f=0;f<nf;++f) {
        if(mesh.facets.nb_vertices(f)!=3) continue;
        k1_max=std::max(k1_max,std::abs(fcurv[f].kappa1));
        k2_max=std::max(k2_max,std::abs(fcurv[f].kappa2));
    }
    auto curv_dist=[&](GEO::index_t a,GEO::index_t b){
        double dk1=(fcurv[a].kappa1-fcurv[b].kappa1)/k1_max;
        double dk2=(fcurv[a].kappa2-fcurv[b].kappa2)/k2_max;
        return dk1*dk1+dk2*dk2;
    };

    std::vector<GEO::index_t> seeds;
    seeds.reserve(actual_count);
    std::vector<double> min_dist(nf,1e30);

    // First seed: most extreme Gaussian curvature
    { double best=-1e30; GEO::index_t s=0;
      for(GEO::index_t f=0;f<nf;++f) {
          if(mesh.facets.nb_vertices(f)!=3) continue;
          double v=std::abs(fcurv[f].gaussian);
          if(v>best){best=v;s=f;}
      }
      seeds.push_back(s); }

    auto update_dists=[&](GEO::index_t seed){
        for(GEO::index_t f=0;f<nf;++f){
            if(mesh.facets.nb_vertices(f)!=3) continue;
            double d=curv_dist(f,seed);
            if(d<min_dist[f]) min_dist[f]=d;
        }
    };
    update_dists(seeds[0]);

    while((GEO::index_t)seeds.size()<actual_count){
        double best=-1.0; GEO::index_t next=seeds[0];
        for(GEO::index_t f=0;f<nf;++f){
            if(mesh.facets.nb_vertices(f)!=3) continue;
            if(min_dist[f]>best){best=min_dist[f];next=f;}
        }
        seeds.push_back(next);
        update_dists(next);
    }

    // BFS growing
    std::queue<GEO::index_t> bfs;
    for(GEO::index_t s=0;s<(GEO::index_t)seeds.size();++s){
        patch_id[seeds[s]]=s; bfs.push(seeds[s]);
    }
    while(!bfs.empty()){
        GEO::index_t f=bfs.front(); bfs.pop();
        GEO::index_t pid=patch_id[f];
        for(GEO::index_t nb:fadj[f])
            if(patch_id[nb]==GEO::index_t(-1)){patch_id[nb]=pid;bfs.push(nb);}
    }
    for(GEO::index_t f=0;f<nf;++f)
        if(patch_id[f]==GEO::index_t(-1)) patch_id[f]=0;

    return actual_count;
}

// =======================================================================
// STAGE 3 — Patch descriptors
// =======================================================================

bool fit_sphere_to_patch(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    GEO::index_t pid,
    PatchDescriptor& desc)
{
    std::unordered_set<GEO::index_t> vset;
    for(GEO::index_t f=0;f<mesh.facets.nb();++f){
        if(patch_id[f]!=pid) continue;
        if(mesh.facets.nb_vertices(f)!=3) continue;
        for(int k=0;k<3;++k) vset.insert(mesh.facets.vertex(f,k));
    }
    std::vector<GEO::vec3> pts;
    pts.reserve(vset.size());
    for(GEO::index_t v:vset) pts.push_back(vertex_point(mesh,v));

    if(pts.size()<4){
        GEO::vec3 cen(0,0,0);
        for(auto& p:pts) cen=cen+p;
        if(!pts.empty()) cen=cen/(double)pts.size();
        desc.sphere_center=cen; desc.sphere_radius=1e6;
        desc.sphere_residual=0.0; return false;
    }
    return fit_sphere_ls(pts,desc.sphere_center,
                         desc.sphere_radius,desc.sphere_residual);
}

void compute_patch_descriptors(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t nb_patches,
    std::vector<PatchDescriptor>& descs)
{
    descs.assign(nb_patches,PatchDescriptor{});
    for(GEO::index_t p=0;p<nb_patches;++p) descs[p].patch_id=p;

    std::vector<double> sg(nb_patches,0),sk1(nb_patches,0),
                        sk2(nb_patches,0),sa(nb_patches,0),
                        se(nb_patches,0);
    std::vector<int> ec(nb_patches,0);

    for(GEO::index_t f=0;f<mesh.facets.nb();++f){
        if(mesh.facets.nb_vertices(f)!=3) continue;
        GEO::index_t pid=patch_id[f];
        if(pid>=nb_patches) continue;
        ++descs[pid].nb_faces;
        sg[pid]+=fcurv[f].gaussian; sk1[pid]+=fcurv[f].kappa1;
        sk2[pid]+=fcurv[f].kappa2;
        GEO::vec3 p0=vertex_point(mesh,mesh.facets.vertex(f,0));
        GEO::vec3 p1=vertex_point(mesh,mesh.facets.vertex(f,1));
        GEO::vec3 p2=vertex_point(mesh,mesh.facets.vertex(f,2));
        sa[pid]+=triangle_area_3pts(p0,p1,p2);
        se[pid]+=GEO::length(p1-p0)+GEO::length(p2-p1)+GEO::length(p0-p2);
        ec[pid]+=3;
    }
    for(GEO::index_t p=0;p<nb_patches;++p){
        int n=descs[p].nb_faces; if(n==0) continue;
        descs[p].mean_gaussian=sg[p]/n; descs[p].mean_kappa1=sk1[p]/n;
        descs[p].mean_kappa2=sk2[p]/n; descs[p].total_area=sa[p];
        descs[p].mean_edge_len=ec[p]>0?se[p]/ec[p]:0.0;
        fit_sphere_to_patch(mesh,patch_id,p,descs[p]);
    }
}

// =======================================================================
// STAGE 4 — K-means
// =======================================================================

double patch_descriptor_distance(
    const PatchDescriptor& a, const PatchDescriptor& b,
    const double inv_std[4])
{
    double d0=(a.sphere_radius -b.sphere_radius )*inv_std[0];
    double d1=(a.mean_kappa1  -b.mean_kappa1  )*inv_std[1];
    double d2=(a.mean_kappa2  -b.mean_kappa2  )*inv_std[2];
    double d3=(a.mean_edge_len-b.mean_edge_len)*inv_std[3];
    return d0*d0+d1*d1+d2*d2+d3*d3;
}

double cluster_patches_kmeans(
    std::vector<PatchDescriptor>& descs,
    GEO::index_t K, int max_iter)
{
    int n=(int)descs.size();
    if(n==0||K==0) return 0.0;
    K=std::min(K,(GEO::index_t)n);

    auto get_feat=[](const PatchDescriptor& d,int i)->double{
        switch(i){case 0:return d.sphere_radius;case 1:return d.mean_kappa1;
                  case 2:return d.mean_kappa2;case 3:return d.mean_edge_len;}
        return 0.0;
    };

    double mean_f[4]={},std_f[4]={};
    for(int fi=0;fi<4;++fi){
        for(auto& d:descs) mean_f[fi]+=get_feat(d,fi);
        mean_f[fi]/=n;
        for(auto& d:descs){double diff=get_feat(d,fi)-mean_f[fi];std_f[fi]+=diff*diff;}
        std_f[fi]=std::sqrt(std_f[fi]/n);
        if(std_f[fi]<1e-15) std_f[fi]=1.0;
    }
    double inv_std[4];
    for(int fi=0;fi<4;++fi) inv_std[fi]=1.0/std_f[fi];

    // K-means++ init
    std::mt19937 rng(42);
    std::vector<int> centre_ids;
    centre_ids.push_back((int)(rng()%(unsigned)n));
    std::vector<double> dists(n,1e30);

    for(GEO::index_t ki=1;ki<K;++ki){
        int last=centre_ids.back();
        double sum=0.0;
        for(int i=0;i<n;++i){
            double d=patch_descriptor_distance(descs[i],descs[last],inv_std);
            if(d<dists[i]) dists[i]=d;
            sum+=dists[i];
        }
        double r=std::uniform_real_distribution<>(0.0,sum)(rng);
        double acc=0.0; int chosen=0;
        for(int i=0;i<n;++i){acc+=dists[i];if(acc>=r){chosen=i;break;}}
        centre_ids.push_back(chosen);
    }

    std::vector<PatchDescriptor> centres(K);
    for(GEO::index_t ki=0;ki<K;++ki) centres[ki]=descs[centre_ids[ki]];

    std::vector<GEO::index_t> assign(n,0);
    double inertia=0.0;

    for(int iter=0;iter<max_iter;++iter){
        inertia=0.0; bool changed=false;
        for(int i=0;i<n;++i){
            double best=1e30; GEO::index_t best_k=0;
            for(GEO::index_t ki=0;ki<K;++ki){
                double d=patch_descriptor_distance(descs[i],centres[ki],inv_std);
                if(d<best){best=d;best_k=ki;}
            }
            if(assign[i]!=best_k) changed=true;
            assign[i]=best_k; inertia+=best;
        }
        if(!changed) break;

        std::vector<double> cnt(K,0),sr(K,0),sk1(K,0),sk2(K,0),se(K,0);
        for(int i=0;i<n;++i){
            GEO::index_t ki=assign[i]; cnt[ki]+=1;
            sr[ki]+=descs[i].sphere_radius; sk1[ki]+=descs[i].mean_kappa1;
            sk2[ki]+=descs[i].mean_kappa2;  se[ki]+=descs[i].mean_edge_len;
        }
        for(GEO::index_t ki=0;ki<K;++ki){
            if(cnt[ki]<1.0) continue;
            centres[ki].sphere_radius=sr[ki]/cnt[ki];
            centres[ki].mean_kappa1=sk1[ki]/cnt[ki];
            centres[ki].mean_kappa2=sk2[ki]/cnt[ki];
            centres[ki].mean_edge_len=se[ki]/cnt[ki];
        }
    }

    for(int i=0;i<n;++i) descs[i].template_id=assign[i];
    return inertia;
}

} // namespace CurveKIDAlgo
} // namespace OGF