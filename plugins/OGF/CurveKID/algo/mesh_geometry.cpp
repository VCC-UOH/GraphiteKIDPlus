// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_geometry.h>
// CVT-based isotropic remeshing (replaces vertex clustering decimator)
#include <geogram/voronoi/CVT.h>
#include <geogram/basic/attributes.h>

#include <OGF/CurveKID/algo/mesh_geometry.h>

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <set>
#include <functional>
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
// STAGE 4 — K-means  (UPDATED: radius cap + log-transform)
// =======================================================================

double patch_descriptor_distance(
    const PatchDescriptor& a, const PatchDescriptor& b,
    const double inv_std[4])
{
    // Feature 0 uses log(r+1) — must match get_feat() in cluster_patches_kmeans
    double d0=(std::log1p(a.sphere_radius) -std::log1p(b.sphere_radius) )*inv_std[0];
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

    // --- Cap degenerate sphere radii before clustering ---
    // Patches with <4 vertices fall back to radius=1e6 (effectively flat).
    // Find the 90th-percentile radius among valid patches and cap at 2x that,
    // so a few tiny/degenerate patches don't collapse the entire feature space.
    {
        std::vector<double> valid_radii;
        valid_radii.reserve(n);
        for (auto& d : descs)
            if (d.nb_faces > 0 && d.sphere_radius < 1e5)
                valid_radii.push_back(d.sphere_radius);

        double cap = 1e5;
        if (!valid_radii.empty()) {
            std::sort(valid_radii.begin(), valid_radii.end());
            int idx = (int)(valid_radii.size() * 0.90);
            if (idx >= (int)valid_radii.size())
                idx = (int)valid_radii.size() - 1;
            cap = valid_radii[idx] * 2.0;
        }
        for (auto& d : descs)
            if (d.sphere_radius > cap) d.sphere_radius = cap;
    }

    // Feature 0: log(sphere_radius + 1) — log-transform squashes the huge
    //            spread between highly-curved (small r) and near-flat (large r)
    //            patches so K-means can distinguish curvature levels properly.
    // Features 1-3: mean_kappa1, mean_kappa2, mean_edge_len (unchanged)
    auto get_feat=[](const PatchDescriptor& d, int i)->double{
        switch(i){
            case 0: return std::log1p(d.sphere_radius);
            case 1: return d.mean_kappa1;
            case 2: return d.mean_kappa2;
            case 3: return d.mean_edge_len;
        }
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

        // Update centroids — recompute as mean of assigned descriptors.
        // Centroid stores raw sphere_radius; get_feat applies log at query time.
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

// =======================================================================
// POTTMANN-INSPIRED INIT — Stage 1: Central spheres
// =======================================================================

void compute_central_spheres(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<CentralSphere>&         cspheres,
    double H_eps)
{
    const GEO::index_t nv = mesh.vertices.nb();
    cspheres.resize(nv);

    // Build 1-ring neighbours for the fallback circumradius
    std::vector<std::vector<GEO::index_t>> ring(nv);
    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
        GEO::index_t nvc = mesh.facets.nb_vertices(f);
        for (GEO::index_t lv = 0; lv < nvc; ++lv) {
            GEO::index_t v  = mesh.facets.vertex(f, lv);
            GEO::index_t v1 = mesh.facets.vertex(f, (lv + 1) % nvc);
            GEO::index_t v2 = mesh.facets.vertex(f, (lv + 2) % nvc);
            ring[v].push_back(v1);
            ring[v].push_back(v2);
        }
    }

    // Pre-compute vertex normals (area-weighted average of face normals)
    std::vector<GEO::vec3> vnormals(nv, GEO::vec3(0, 0, 0));
    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
        if (mesh.facets.nb_vertices(f) < 3) continue;
        GEO::vec3 p0 = vertex_point(mesh, mesh.facets.vertex(f, 0));
        GEO::vec3 p1 = vertex_point(mesh, mesh.facets.vertex(f, 1));
        GEO::vec3 p2 = vertex_point(mesh, mesh.facets.vertex(f, 2));
        GEO::vec3 fn = GEO::cross(p1 - p0, p2 - p0); // magnitude = 2*area
        for (GEO::index_t lv = 0; lv < mesh.facets.nb_vertices(f); ++lv)
            vnormals[mesh.facets.vertex(f, lv)] = vnormals[mesh.facets.vertex(f, lv)] + fn;
    }
    for (GEO::index_t v = 0; v < nv; ++v) {
        double len = GEO::length(vnormals[v]);
        if (len > 1e-15) vnormals[v] = vnormals[v] * (1.0 / len);
    }

    for (GEO::index_t v = 0; v < nv; ++v) {
        CentralSphere& cs = cspheres[v];
        cs.H_value = vcurv[v].mean;

        double absH = std::fabs(cs.H_value);
        if (absH >= H_eps) {
            // Normal case: centre = vertex + (1/H) * normal
            double inv_H = 1.0 / cs.H_value;
            GEO::vec3 vp = vertex_point(mesh, v);
            cs.center = vp + vnormals[v] * inv_H;
            cs.radius = std::fabs(inv_H);
            cs.valid  = true;
        } else {
            // Near-flat or saddle near H=0:
            // Fallback: estimate radius as mean distance to 1-ring neighbours
            // (proxy for local curvature radius when H≈0)
            GEO::vec3 vp = vertex_point(mesh, v);
            double sum_dist = 0.0;
            int    cnt      = 0;
            // deduplicate ring
            std::unordered_set<GEO::index_t> seen;
            for (GEO::index_t nb : ring[v]) {
                if (!seen.insert(nb).second) continue;
                sum_dist += GEO::length(vertex_point(mesh, nb) - vp);
                ++cnt;
            }
            double fallback_r = (cnt > 0) ? (sum_dist / cnt) : 1.0;
            // Centre: vertex + fallback_r * normal (approximates flat sphere)
            cs.center = vp + vnormals[v] * fallback_r;
            cs.radius = fallback_r;
            cs.valid  = false;  // flagged so callers can filter
        }
    }
}

// =======================================================================
// POTTMANN-INSPIRED INIT — Stage 2: Curvature-adaptive target edge lengths
// =======================================================================

void compute_target_edge_lengths(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<double>&                target_len,
    double H_eps)
{
    const GEO::index_t nf = mesh.facets.nb();
    target_len.resize(nf, 0.0);

    // g(f) = c / sqrt( max(|H_f|, H_eps) )
    // c chosen so mean g = mesh mean edge length
    // Pass 1: collect raw g values
    std::vector<double> raw(nf);
    for (GEO::index_t f = 0; f < nf; ++f) {
        double H_sum = 0.0;
        int cnt = 0;
        for (GEO::index_t lv = 0; lv < mesh.facets.nb_vertices(f); ++lv) {
            H_sum += std::fabs(vcurv[mesh.facets.vertex(f, lv)].mean);
            ++cnt;
        }
        double H_face = (cnt > 0) ? H_sum / cnt : H_eps;
        raw[f] = 1.0 / std::sqrt(std::max(H_face, H_eps));
    }

    // Compute mesh mean edge length for normalisation
    double total_edge = 0.0;
    int    edge_cnt   = 0;
    for (GEO::index_t f = 0; f < nf; ++f) {
        GEO::index_t nvc = mesh.facets.nb_vertices(f);
        for (GEO::index_t lv = 0; lv < nvc; ++lv) {
            GEO::vec3 a = vertex_point(mesh, mesh.facets.vertex(f, lv));
            GEO::vec3 b = vertex_point(mesh, mesh.facets.vertex(f, (lv+1) % nvc));
            total_edge += GEO::length(b - a);
            ++edge_cnt;
        }
    }
    double mean_edge = (edge_cnt > 0) ? total_edge / edge_cnt : 1.0;

    // Pass 2: scale so mean(g) = mean_edge
    double mean_raw = 0.0;
    for (double r : raw) mean_raw += r;
    mean_raw /= (nf > 0 ? nf : 1);
    double c = (mean_raw > 1e-15) ? mean_edge / mean_raw : 1.0;

    // Clamp to [0.25*mean, 4*mean]
    double lo = 0.25 * mean_edge;
    double hi = 4.0  * mean_edge;
    for (GEO::index_t f = 0; f < nf; ++f)
        target_len[f] = std::max(lo, std::min(hi, c * raw[f]));
}

// =======================================================================
// POTTMANN-INSPIRED INIT — Stage 3: Fit one SphericalFace per triangle
// =======================================================================

void initialize_spherical_faces(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    const std::vector<CentralSphere>&   cspheres,
    std::vector<SphericalFace>&         sfaces)
{
    const GEO::index_t nf = mesh.facets.nb();
    sfaces.resize(nf);

    for (GEO::index_t f = 0; f < nf; ++f) {
        SphericalFace& sf = sfaces[f];
        if (mesh.facets.nb_vertices(f) < 3) continue;

        sf.v[0] = mesh.facets.vertex(f, 0);
        sf.v[1] = mesh.facets.vertex(f, 1);
        sf.v[2] = mesh.facets.vertex(f, 2);

        GEO::vec3 p[3];
        for (int k = 0; k < 3; ++k) p[k] = vertex_point(mesh, sf.v[k]);

        // --- Fit sphere: least-squares to 3 vertex positions + their normals
        // as tangent-plane constraints.
        //
        // Three vertices alone don't determine a sphere uniquely (underdetermined).
        // We add the constraint: the sphere-centre-to-vertex vector is parallel
        // to the vertex normal, i.e. (C - p_i) × n_i = 0 for each i.
        //
        // Practical approach: pick the vertex with largest |H| as the "primary"
        // central sphere (from cspheres), then refine the centre by pulling
        // it toward the centroid of the three corresponding central sphere centres.

        // Identify the most-curved corner (best-conditioned central sphere)
        int best_k = 0;
        double best_H = std::fabs(vcurv[sf.v[0]].mean);
        for (int k = 1; k < 3; ++k) {
            double h = std::fabs(vcurv[sf.v[k]].mean);
            if (h > best_H) { best_H = h; best_k = k; }
        }

        // Weighted average of the three central sphere centres
        // weight = |H_i|  (more-curved corners get more influence)
        GEO::vec3 C_avg(0, 0, 0);
        double w_sum = 0.0;
        for (int k = 0; k < 3; ++k) {
            double w = std::fabs(vcurv[sf.v[k]].mean);
            if (!cspheres[sf.v[k]].valid) w *= 0.1; // downweight flat-fallback
            C_avg = C_avg + cspheres[sf.v[k]].center * w;
            w_sum += w;
        }
        if (w_sum < 1e-15) {
            // All three vertices near-flat: use circumcenter of triangle,
            // pushed along face normal to match best corner radius
            GEO::vec3 edge1 = p[1] - p[0];
            GEO::vec3 edge2 = p[2] - p[0];
            GEO::vec3 fn    = GEO::cross(edge1, edge2);
            double fn_len = GEO::length(fn);
            GEO::vec3 centroid = (p[0] + p[1] + p[2]) * (1.0 / 3.0);
            double fallback_r = cspheres[sf.v[best_k]].radius;
            GEO::vec3 fn_unit = (fn_len > 1e-15) ? fn * (1.0 / fn_len) : GEO::vec3(0,0,1);
            C_avg = centroid + fn_unit * fallback_r;
            w_sum = 1.0;
        }
        GEO::vec3 C = C_avg * (1.0 / w_sum);

        // Radius = mean distance from C to 3 vertices
        double R = 0.0;
        for (int k = 0; k < 3; ++k) R += GEO::length(p[k] - C);
        R /= 3.0;

        // RMS residual
        double rms = 0.0;
        for (int k = 0; k < 3; ++k) {
            double d = GEO::length(p[k] - C) - R;
            rms += d * d;
        }
        sf.sphere_center = C;
        sf.sphere_radius  = R;
        sf.fit_residual   = std::sqrt(rms / 3.0);
        sf.valid          = (R > 1e-10);

        // --- Pose-invariant descriptor
        sf.desc_radius = R;

        // Arc lengths: great-circle distance = R * central_angle
        for (int k = 0; k < 3; ++k) {
            GEO::vec3 a = p[k]        - C;
            GEO::vec3 b = p[(k+1)%3]  - C;
            double la = GEO::length(a);
            double lb = GEO::length(b);
            double arc = 0.0;
            if (la > 1e-15 && lb > 1e-15) {
                double cos_theta = GEO::dot(a, b) / (la * lb);
                cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
                arc = R * std::acos(cos_theta);
            }
            sf.desc_arc[k] = arc;
        }

        // Solid angle subtended by triangle at sphere centre
        {
            GEO::vec3 a = p[0] - C, b = p[1] - C, c = p[2] - C;
            double la = GEO::length(a), lb = GEO::length(b), lc = GEO::length(c);
            double solid_angle = 0.0;
            if (la > 1e-15 && lb > 1e-15 && lc > 1e-15) {
                // Van Oosterom & Strackee formula
                double num = std::fabs(GEO::dot(a, GEO::cross(b, c)));
                double denom = la*lb*lc
                    + GEO::dot(a,b)*lc
                    + GEO::dot(b,c)*la
                    + GEO::dot(c,a)*lb;
                solid_angle = 2.0 * std::atan2(num, denom);
                if (solid_angle < 0.0) solid_angle += 4.0 * M_PI; // wrap
            }
            sf.desc_solid_angle = solid_angle;
        }

        // Mean |H| over 3 vertices
        sf.desc_mean_H = (std::fabs(vcurv[sf.v[0]].mean)
                        + std::fabs(vcurv[sf.v[1]].mean)
                        + std::fabs(vcurv[sf.v[2]].mean)) / 3.0;
    }
}

// =======================================================================
// CURVED MESH CONSTRUCTION — Stage 1: coarse remesh
// =======================================================================

// =======================================================================
// TOPOLOGY VALIDATION
// =======================================================================

MeshTopologyReport validate_mesh_topology(const GEO::Mesh& mesh, bool log_results)
{
    MeshTopologyReport r;
    const GEO::index_t nv = mesh.vertices.nb();
    const GEO::index_t nf = mesh.facets.nb();

    // Count edges and their valence using a map {min_v, max_v} -> count
    std::map<std::pair<GEO::index_t,GEO::index_t>, int> edge_valence;
    int degenerate = 0;
    for (GEO::index_t f = 0; f < nf; ++f) {
        if (mesh.facets.nb_vertices(f) != 3) continue;
        GEO::index_t v[3];
        v[0] = mesh.facets.vertex(f, 0);
        v[1] = mesh.facets.vertex(f, 1);
        v[2] = mesh.facets.vertex(f, 2);
        if (v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) { ++degenerate; continue; }
        for (int e = 0; e < 3; ++e) {
            GEO::index_t a = v[e], b = v[(e+1)%3];
            if (a > b) std::swap(a, b);
            edge_valence[std::make_pair(a, b)]++;
        }
    }
    r.degenerate_faces = degenerate;

    int ne = (int)edge_valence.size();
    int boundary_e = 0, nonmanifold_e = 0;
    for (auto& kv : edge_valence) {
        if (kv.second == 1) ++boundary_e;
        if (kv.second > 2) ++nonmanifold_e;
    }
    r.boundary_edges    = boundary_e;
    r.nonmanifold_edges = nonmanifold_e;

    // Euler characteristic
    r.euler = (int)nv - ne + (int)nf;

    // Connected components via union-find on vertices
    std::vector<GEO::index_t> par(nv);
    std::iota(par.begin(), par.end(), 0);
    std::function<GEO::index_t(GEO::index_t)> root = [&](GEO::index_t x) -> GEO::index_t {
        while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x;
    };
    for (GEO::index_t f = 0; f < nf; ++f) {
        if (mesh.facets.nb_vertices(f) != 3) continue;
        GEO::index_t ra = root(mesh.facets.vertex(f,0));
        GEO::index_t rb = root(mesh.facets.vertex(f,1));
        GEO::index_t rc = root(mesh.facets.vertex(f,2));
        par[ra] = rb; par[root(rb)] = root(rc);
    }
    std::set<GEO::index_t> roots;
    for (GEO::index_t v = 0; v < nv; ++v) roots.insert(root(v));
    r.connected_components = (int)roots.size();

    r.face_count_even = (nf % 2 == 0);
    // Accept any connected manifold — closed or open.
    // Only reject non-manifold meshes and multiple disconnected components.
    // FaceCountEven is informational only — not required for open surfaces.
    // Accept any manifold mesh — single or multiple components, open or closed.
    // The sphere fitting works per-face regardless of global topology.
    r.valid = (r.nonmanifold_edges == 0 &&
               r.degenerate_faces == 0);

    if (log_results) {
        GEO::Logger::out("CurveKID") << "  --- Topology validation ---" << std::endl;
        GEO::Logger::out("CurveKID") << "  V=" << nv << " E=" << ne << " F=" << nf
            << "  Euler=" << r.euler
            << "  BoundaryEdges=" << r.boundary_edges << std::endl;
        GEO::Logger::out("CurveKID") << "  Components=" << r.connected_components
            << "  BoundaryEdges=" << r.boundary_edges
            << "  NonmanifoldEdges=" << r.nonmanifold_edges
            << "  DegenerateFaces=" << r.degenerate_faces << std::endl;
        GEO::Logger::out("CurveKID") << "  FaceCountEven=" << (r.face_count_even?"yes":"NO")
            << "  Valid=" << (r.valid?"YES":"NO") << std::endl;
    }
    return r;
}

// Forward declaration of the internal decimator (used only as fallback)
static GEO::index_t build_coarse_mesh_decimate(
    const GEO::Mesh&                    original,
    const std::vector<VertexCurvature>& vcurv,
    GEO::Mesh&                          coarse_mesh,
    GEO::index_t                        target_faces,
    double                              H_eps);

GEO::index_t build_coarse_mesh(
    const GEO::Mesh&                    original,
    const std::vector<VertexCurvature>& vcurv,
    GEO::Mesh&                          coarse_mesh,
    GEO::index_t                        target_faces,
    double                              H_eps)
{
    // CVT-based isotropic remeshing (replaces vertex clustering decimator).
    // For a closed surface with F faces, we need approximately F/2 vertices,
    // so nb_points = target_faces / 2.
    //
    // GEO::CentroidalVoronoiTesselation:
    //   1. compute_initial_sampling(nb_points) — place points uniformly
    //   2. Lloyd_iterations(nb_iter)           — relax to CVT
    //   3. compute_surface(coarse_mesh)        — build output triangulation
    //
    // This produces near-equilateral triangles with no seam artifacts —
    // exactly what we need for reliable sphere fitting.

    // V ≈ F/2 for genus-0 (sphere): Euler V-E+F=2, E=3F/2 → V=F/2+1.
    // For genus-1 (torus): Euler=0 → V=F/2. Same formula works for both.
    // CVT produces approximately 2*nb_points faces for a closed surface.
    GEO::index_t nb_points = std::max(GEO::index_t(6), target_faces / 2);

    GEO::Logger::out("CurveKID")
        << "  CVT remesh: nb_points=" << nb_points
        << "  target_faces=" << target_faces << std::endl;

    // Work on a copy so we don't modify original
    GEO::Mesh src;
    src.copy(original, false);
    GEO::mesh_repair(src, GEO::MESH_REPAIR_TRIANGULATE);

    GEO::CentroidalVoronoiTesselation cvt(&src);
    cvt.set_use_RVC_centroids(true);  // matches remesh:RVC_centroids=true default
    cvt.compute_initial_sampling(nb_points);
    cvt.Lloyd_iterations(40);         // matches opt:nb_Lloyd_iter=40 (no HLBFGS)
    cvt.compute_surface(&coarse_mesh, true); // multinerve=true preserves topology

    GEO::mesh_repair(coarse_mesh, GEO::MESH_REPAIR_DEFAULT);

    GEO::Logger::out("CurveKID")
        << "  CVT result: " << coarse_mesh.facets.nb() << " faces, "
        << coarse_mesh.vertices.nb() << " vertices" << std::endl;

    return coarse_mesh.facets.nb();
}

// =======================================================================
// CURVED MESH CONSTRUCTION — Stage 1 (fallback): hand-written decimator
// Kept for reference; build_coarse_mesh now uses mesh_decimate_full above.
// =======================================================================
static GEO::index_t build_coarse_mesh_decimate(
    const GEO::Mesh&                    original,
    const std::vector<VertexCurvature>& vcurv,
    GEO::Mesh&                          coarse_mesh,
    GEO::index_t                        target_faces,
    double                              H_eps)
{
    // Compute mean edge length of original mesh
    double total_len = 0.0;
    int    edge_cnt  = 0;
    for (GEO::index_t f = 0; f < original.facets.nb(); ++f) {
        GEO::index_t nvc = original.facets.nb_vertices(f);
        for (GEO::index_t lv = 0; lv < nvc; ++lv) {
            GEO::vec3 a = vertex_point(original, original.facets.vertex(f, lv));
            GEO::vec3 b = vertex_point(original, original.facets.vertex(f, (lv+1)%nvc));
            total_len += GEO::length(b - a);
            ++edge_cnt;
        }
    }
    double mean_edge = (edge_cnt > 0) ? total_len / edge_cnt : 1.0;

    // Target edge length for coarse mesh:
    // A mesh of T equilateral triangles on area A has edge ~ sqrt(2A/sqrt(3)/T)
    // Simpler: scale mean_edge by sqrt(original_faces / target_faces)
    double orig_faces = double(original.facets.nb());
    double scale      = (target_faces > 0 && orig_faces > 0)
                        ? std::sqrt(orig_faces / double(target_faces))
                        : 1.0;
    double target_edge = mean_edge * scale;

    // -----------------------------------------------------------------------
    // Self-contained greedy edge-collapse decimation.
    // No dependency on mesh_remesh_smooth or mesh_decimate_full —
    // both are absent in some Geogram distributions.
    //
    // Strategy: repeatedly collapse the edge whose midpoint moves the
    // surface the least (weighted by curvature so flat regions collapse
    // first).  Stop when face count <= target_faces.
    // -----------------------------------------------------------------------

    // 1. Copy original into coarse_mesh (vertices + triangle connectivity only)
    coarse_mesh.clear();
    {
        const GEO::index_t nv = original.vertices.nb();
        const GEO::index_t nf = original.facets.nb();
        coarse_mesh.vertices.create_vertices(nv);
        for (GEO::index_t v = 0; v < nv; ++v) {
            double* dst = coarse_mesh.vertices.point_ptr(v);
            const double* src = original.vertices.point_ptr(v);
            dst[0] = src[0]; dst[1] = src[1]; dst[2] = src[2];
        }
        for (GEO::index_t f = 0; f < nf; ++f) {
            if (original.facets.nb_vertices(f) == 3)
                coarse_mesh.facets.create_triangle(
                    original.facets.vertex(f,0),
                    original.facets.vertex(f,1),
                    original.facets.vertex(f,2));
        }
        coarse_mesh.facets.connect();
    }

    // 2. Compute per-vertex curvature weight for original mesh
    //    (used to bias collapse cost: low-curvature vertices collapse first)
    const GEO::index_t nv_orig = original.vertices.nb();
    std::vector<double> curv_weight(nv_orig, 1.0);
    for (GEO::index_t v = 0; v < nv_orig; ++v) {
        double H = std::fabs(vcurv[v].mean);
        curv_weight[v] = std::sqrt(std::max(H, H_eps) / H_eps);
    }

    // 3. Greedy collapse loop
    // We track which vertices/faces are still active
    const GEO::index_t NV = coarse_mesh.vertices.nb();
    const GEO::index_t NF_orig = coarse_mesh.facets.nb();

    std::vector<bool> v_active(NV, true);
    std::vector<bool> f_active(NF_orig, true);
    // parent[v] = canonical representative after collapse
    std::vector<GEO::index_t> parent(NV);
    std::iota(parent.begin(), parent.end(), 0);

    // Cache face vertex arrays so we never call root() during deactivation
    // after parent[] has already been updated — that was the original bug.
    // fv[f][k] stores the ORIGINAL vertex index; we resolve through root()
    // only when we need the current representative, and we do so before
    // any parent[] update in that iteration.
    // Store face vertices as flat array: fv[f*3+k]
    // Use size_t indexing throughout to avoid MSVC signed/unsigned mix errors.
    std::vector<GEO::index_t> fv((size_t)NF_orig * 3u);
    for (size_t f = 0; f < (size_t)NF_orig; ++f)
        for (size_t k = 0; k < 3u; ++k)
            fv[f*3u+k] = coarse_mesh.facets.vertex((GEO::index_t)f, (GEO::index_t)k);

    // Build vertex → incident faces adjacency for O(deg) neighbour queries
    std::vector<std::vector<GEO::index_t>> v2f((size_t)NV);
    for (size_t f = 0; f < (size_t)NF_orig; ++f)
        for (size_t k = 0; k < 3u; ++k)
            v2f[fv[f*3u+k]].push_back((GEO::index_t)f);

    auto root_fn = [&](GEO::index_t v) -> GEO::index_t {
        while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
        return v;
    };

    // Helper: get current vertex position (goes through root)
    auto vpos = [&](GEO::index_t v) -> GEO::vec3 {
        const double* p = coarse_mesh.vertices.point_ptr(root_fn(v));
        return GEO::vec3(p[0], p[1], p[2]);
    };

    GEO::index_t active_faces = NF_orig;

    // Priority queue: (cost, v0, v1) — collapse edge v0-v1 to midpoint
    using EdgeCost = std::tuple<double, GEO::index_t, GEO::index_t>;
    std::priority_queue<EdgeCost, std::vector<EdgeCost>, std::greater<EdgeCost>> pq;

    // Seed queue with all edges
    {
        std::set<std::pair<GEO::index_t,GEO::index_t>> seen_edges;
        for (size_t f = 0; f < (size_t)NF_orig; ++f) {
            for (size_t lv = 0; lv < 3u; ++lv) {
                GEO::index_t a = fv[f*3u+lv], b = fv[f*3u+(lv+1u)%3u];
                if (a > b) std::swap(a, b);
                if (seen_edges.insert(std::make_pair(a,b)).second) {
                    double len = GEO::length(vpos(b) - vpos(a));
                    double w = (curv_weight[a] + curv_weight[b]) * 0.5;
                    // Multiply by w: high-curvature edges get HIGHER cost → preserved longer
                    pq.push(EdgeCost(len * w, a, b));
                }
            }
        }
    }

    while (active_faces > target_faces && !pq.empty()) {
        EdgeCost top_ec = pq.top(); pq.pop();
        double cost = std::get<0>(top_ec);
        GEO::index_t va = std::get<1>(top_ec), vb = std::get<2>(top_ec);
        (void)cost;

        // Resolve current representatives BEFORE any parent[] update
        GEO::index_t ra = root_fn(va), rb = root_fn(vb);
        if (ra == rb) continue;          // already merged
        if (!v_active[ra] || !v_active[rb]) continue;

        // Move ra to midpoint
        GEO::vec3 pa = vpos(ra), pb = vpos(rb);
        GEO::vec3 mid = (pa + pb) * 0.5;
        {
            double* p = coarse_mesh.vertices.point_ptr(ra);
            p[0] = mid.x; p[1] = mid.y; p[2] = mid.z;
        }

        // Merge rb → ra AFTER we have captured ra/rb
        parent[rb] = ra;
        v_active[rb] = false;

        // Deactivate faces that have become degenerate after the merge.
        // A face is degenerate if any two of its three resolved vertices are equal.
        // We scan ALL active faces (correctness-first) — O(NF) per collapse but
        // guarantees no valid face is wrongly removed.
        // Note: must call root_fn AFTER parent[rb]=ra so rb resolves to ra.
        for (GEO::index_t fi = 0; fi < NF_orig; ++fi) {
            if (!f_active[fi]) continue;
            GEO::index_t va2 = root_fn(fv[(size_t)fi*3u]);
            GEO::index_t vb2 = root_fn(fv[(size_t)fi*3u+1]);
            GEO::index_t vc2 = root_fn(fv[(size_t)fi*3u+2]);
            if (va2 == vb2 || vb2 == vc2 || va2 == vc2) {
                f_active[fi] = false;
                --active_faces;
            }
        }

        // Add new edges from ra to its current neighbours (faces incident to ra
        // after the merge, i.e. union of v2f[ra] and v2f[rb])
        // Enqueue new edges from ra to neighbours (union of ra and rb incident faces)
        for (int pass = 0; pass < 2; ++pass) {
            GEO::index_t v_orig = (pass == 0) ? ra : rb;
            for (size_t fi = 0; fi < v2f[v_orig].size(); ++fi) {
                GEO::index_t f = v2f[v_orig][fi];
                if (!f_active[f]) continue;
                for (size_t k = 0; k < 3u; ++k) {
                    GEO::index_t rn = root_fn(fv[(size_t)f*3u+k]);
                    if (rn == ra) continue;
                    GEO::index_t ea = std::min(ra, rn), eb = std::max(ra, rn);
                    double len = GEO::length(vpos(rn) - vpos(ra));
                    double w = (curv_weight[ra] + curv_weight[rn]) * 0.5;
                    pq.push(EdgeCost(len * w, ea, eb));
                }
            }
        }
    }

    // 4. Rebuild coarse_mesh from active faces / remapped vertices
    {
        // Remap active vertices to compact indices
        std::vector<GEO::index_t> new_idx(NV, GEO::index_t(-1));
        GEO::Mesh result;
        for (GEO::index_t v = 0; v < NV; ++v) {
            if (v_active[v]) {
                new_idx[v] = result.vertices.create_vertices(1);
                const double* src = coarse_mesh.vertices.point_ptr(v);
                double* dst = result.vertices.point_ptr(new_idx[v]);
                dst[0]=src[0]; dst[1]=src[1]; dst[2]=src[2];
            }
        }
        for (size_t f = 0; f < (size_t)NF_orig; ++f) {
            if (!f_active[(GEO::index_t)f]) continue;
            GEO::index_t rv[3];
            for (size_t k = 0; k < 3u; ++k)
                rv[k] = new_idx[root_fn(fv[f*3u+k])];
            if (rv[0] == rv[1] || rv[1] == rv[2] || rv[0] == rv[2]) continue;
            if (rv[0]==GEO::index_t(-1)||rv[1]==GEO::index_t(-1)||rv[2]==GEO::index_t(-1)) continue;
            result.facets.create_triangle(rv[0], rv[1], rv[2]);
        }
        result.facets.connect();
        coarse_mesh.copy(result, false);
    }

    return coarse_mesh.facets.nb();
}

// =======================================================================
// CURVED MESH CONSTRUCTION — Stage 2: fit sphere per coarse triangle
// =======================================================================

void fit_curved_faces(
    const GEO::Mesh&         original,
    const GEO::Mesh&         coarse_mesh,
    std::vector<CurvedFace>& cfaces)
{
    const GEO::index_t nf_coarse = coarse_mesh.facets.nb();
    cfaces.resize(nf_coarse);

    // --- Step 0: compute vertex normals + mean curvature on the ORIGINAL mesh ---
    const GEO::index_t nv_orig = original.vertices.nb();
    std::vector<GEO::vec3> orig_pts(nv_orig);
    std::vector<GEO::vec3> orig_normals(nv_orig, GEO::vec3(0,0,0));
    std::vector<double>    orig_H(nv_orig, 0.0);  // mean curvature magnitude

    for (GEO::index_t v = 0; v < nv_orig; ++v)
        orig_pts[v] = vertex_point(original, v);

    // Area-weighted vertex normals
    for (GEO::index_t f = 0; f < original.facets.nb(); ++f) {
        if (original.facets.nb_vertices(f) < 3) continue;
        GEO::vec3 a = orig_pts[original.facets.vertex(f,0)];
        GEO::vec3 b = orig_pts[original.facets.vertex(f,1)];
        GEO::vec3 c = orig_pts[original.facets.vertex(f,2)];
        GEO::vec3 fn = GEO::cross(b-a, c-a);
        for (GEO::index_t lv = 0; lv < 3; ++lv)
            orig_normals[original.facets.vertex(f,lv)] =
                orig_normals[original.facets.vertex(f,lv)] + fn;
    }
    for (GEO::index_t v = 0; v < nv_orig; ++v) {
        double len = GEO::length(orig_normals[v]);
        if (len > 1e-14) orig_normals[v] = orig_normals[v] * (1.0/len);
    }

    // Per-vertex mean curvature |H| via cotangent formula
    // H = |laplacian(x)| / 2  (Laplace-Beltrami of position)
    std::vector<double> vert_area(nv_orig, 0.0);
    std::vector<GEO::vec3> laplacian(nv_orig, GEO::vec3(0,0,0));
    for (GEO::index_t f = 0; f < original.facets.nb(); ++f) {
        if (original.facets.nb_vertices(f) < 3) continue;
        GEO::index_t vi[3];
        GEO::vec3    vp[3];
        for (int k = 0; k < 3; ++k) {
            vi[k] = original.facets.vertex(f, k);
            vp[k] = orig_pts[vi[k]];
        }
        double face_area = 0.5 * GEO::length(GEO::cross(vp[1]-vp[0], vp[2]-vp[0]));
        for (int k = 0; k < 3; ++k) {
            vert_area[vi[k]] += face_area / 3.0;
            // cotangent weight for opposite edge
            int i = k, j = (k+1)%3, m = (k+2)%3;
            GEO::vec3 e1 = vp[i] - vp[m], e2 = vp[j] - vp[m];
            double cross_len = GEO::length(GEO::cross(e1, e2));
            double cot = (cross_len > 1e-14)
                       ? GEO::dot(e1, e2) / cross_len : 0.0;
            laplacian[vi[i]] = laplacian[vi[i]] + (vp[j] - vp[i]) * cot;
            laplacian[vi[j]] = laplacian[vi[j]] + (vp[i] - vp[j]) * cot;
        }
    }
    for (GEO::index_t v = 0; v < nv_orig; ++v) {
        if (vert_area[v] > 1e-14)
            orig_H[v] = 0.5 * GEO::length(laplacian[v]) / vert_area[v];
    }

    // Global mean and std of |H| — used to set curvature threshold adaptively
    double H_mean = 0.0;
    for (GEO::index_t v = 0; v < nv_orig; ++v) H_mean += orig_H[v];
    H_mean /= std::max((GEO::index_t)1, nv_orig);

    // H_thresh unused — closest-point approach never rejects by curvature
    (void)H_mean;  // suppress unused warning

    // --- Step 1a: classify each coarse face as smooth or seam ---
    // A face is a "seam face" if any of its edges has a dihedral angle > SHARP_DEG
    // with the adjacent face. Seam faces straddle sharp features and produce
    // corrupted sphere fits — they are kept flat (valid=false).
    // We also dilate: immediate neighbours of seam faces are also kept flat,
    // to avoid sphere fits that collect points from across the seam.
    const double SHARP_DEG = 175.0;  // effectively disabled — seam detection off
    const double SHARP_COS = std::cos(SHARP_DEG * 3.14159265358979323846 / 180.0);

    // Compute coarse face normals
    std::vector<GEO::vec3> coarse_fnormals(nf_coarse, GEO::vec3(0,0,0));
    for (GEO::index_t f = 0; f < nf_coarse; ++f) {
        if (coarse_mesh.facets.nb_vertices(f) < 3) continue;
        GEO::vec3 a = vertex_point(coarse_mesh, coarse_mesh.facets.vertex(f,0));
        GEO::vec3 b = vertex_point(coarse_mesh, coarse_mesh.facets.vertex(f,1));
        GEO::vec3 c = vertex_point(coarse_mesh, coarse_mesh.facets.vertex(f,2));
        GEO::vec3 n = GEO::cross(b-a, c-a);
        double len = GEO::length(n);
        if (len > 1e-14) n = n * (1.0/len);
        coarse_fnormals[f] = n;
    }

    // Mark seam faces: any adjacent edge with dihedral angle > SHARP_DEG
    std::vector<bool> is_seam(nf_coarse, false);
    for (GEO::index_t f = 0; f < nf_coarse; ++f) {
        GEO::index_t nle = coarse_mesh.facets.nb_vertices(f);
        for (GEO::index_t le = 0; le < nle; ++le) {
            GEO::index_t adj = coarse_mesh.facets.adjacent(f, le);
            if (adj == GEO::NO_FACET || adj >= nf_coarse) continue;
            double d = GEO::dot(coarse_fnormals[f], coarse_fnormals[adj]);
            // d < SHARP_COS means angle > SHARP_DEG
            if (d < SHARP_COS) { is_seam[f] = true; break; }
        }
    }

    // No dilation — only the seam faces themselves are kept flat.
    // For smooth surfaces there are no true seams, so dilation would
    // incorrectly mark valid smooth faces as flat.

    // Count seam vs smooth faces for logging
    int seam_cnt = 0;
    for (GEO::index_t f = 0; f < nf_coarse; ++f) if (is_seam[f]) ++seam_cnt;
    GEO::Logger::out("CurveKID")
        << "  Sharp-feature classification: " << seam_cnt << " seam/flat faces, "
        << (nf_coarse - seam_cnt) << " smooth faces (threshold " << SHARP_DEG << " deg)"
        << std::endl;

    // --- Step 1: fit one sphere per coarse face ---
    for (GEO::index_t f = 0; f < nf_coarse; ++f) {
        CurvedFace& cf = cfaces[f];
        cf.face_id = f;
        cf.valid   = false;

        if (coarse_mesh.facets.nb_vertices(f) < 3) continue;

        // Skip seam faces
        if (is_seam[f]) continue;

        // Skip boundary faces — their normals are unreliable (no neighbour on one side)
        // These faces produce spikes at open mesh boundaries.
        bool is_boundary_face = false;
        for (GEO::index_t le = 0; le < coarse_mesh.facets.nb_vertices(f); ++le) {
            if (coarse_mesh.facets.adjacent(f, le) == GEO::NO_FACET) {
                is_boundary_face = true; break;
            }
        }
        if (is_boundary_face) continue;

        cf.v[0] = coarse_mesh.facets.vertex(f, 0);
        cf.v[1] = coarse_mesh.facets.vertex(f, 1);
        cf.v[2] = coarse_mesh.facets.vertex(f, 2);

        GEO::vec3 p[3];
        for (int k = 0; k < 3; ++k)
            p[k] = vertex_point(coarse_mesh, cf.v[k]);

        GEO::vec3 centroid = (p[0]+p[1]+p[2])*(1.0/3.0);

        // Circumradius of this coarse triangle
        double ea = GEO::length(p[1]-p[0]);
        double eb = GEO::length(p[2]-p[1]);
        double ec = GEO::length(p[0]-p[2]);
        double hs = (ea+eb+ec)*0.5;
        double tri_area = std::sqrt(std::max(0.0, hs*(hs-ea)*(hs-eb)*(hs-ec)));
        double circum_r = (tri_area > 1e-14)
                        ? (ea*eb*ec)/(4.0*tri_area)
                        : (ea+eb+ec)/3.0;

        // Coarse face normal
        GEO::vec3 coarse_normal = GEO::cross(p[1]-p[0], p[2]-p[0]);
        double cn_len = GEO::length(coarse_normal);
        if (cn_len > 1e-14) coarse_normal = coarse_normal * (1.0/cn_len);

        // ---------------------------------------------------------------
        // POTTMANN Eq.7 — CORRECT IMPLEMENTATION
        //
        // For each coarse vertex, find its closest original vertex and
        // use that original vertex's position + normal in Eq.7.
        // This is exactly what Pottmann does: normals at the face corners.
        // No neighbourhood search, no radius parameter, no sliver rejection.
        // ---------------------------------------------------------------
        std::vector<GEO::vec3> nbr_pts;
        std::vector<GEO::vec3> nbr_nrm;

        // For each coarse corner: find closest original vertex
        for (int k = 0; k < 3; ++k) {
            double best_d2 = std::numeric_limits<double>::max();
            GEO::index_t best_v = 0;
            for (GEO::index_t v = 0; v < nv_orig; ++v) {
                GEO::vec3 diff = orig_pts[v] - p[k];
                double d2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
                if (d2 < best_d2) { best_d2 = d2; best_v = v; }
            }
            nbr_pts.push_back(orig_pts[best_v]);
            nbr_nrm.push_back(orig_normals[best_v]);
        }

        // Also sample interior: face centroid and edge midpoints → 4 more points
        // This over-constrains the 3×3 system for a more stable solution.
        GEO::vec3 sample_pts[4] = {
            centroid,
            (p[0]+p[1])*0.5,
            (p[1]+p[2])*0.5,
            (p[2]+p[0])*0.5
        };
        for (int k = 0; k < 4; ++k) {
            double best_d2 = std::numeric_limits<double>::max();
            GEO::index_t best_v = 0;
            for (GEO::index_t v = 0; v < nv_orig; ++v) {
                GEO::vec3 diff = orig_pts[v] - sample_pts[k];
                double d2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
                if (d2 < best_d2) { best_d2 = d2; best_v = v; }
            }
            nbr_pts.push_back(orig_pts[best_v]);
            nbr_nrm.push_back(orig_normals[best_v]);
        }

        // Compute mean H for descriptor (informational only — not used for rejection)
        double mean_H_local = 0.0;
        for (GEO::index_t v = 0; v < nv_orig; ++v) {
            GEO::vec3 diff = orig_pts[v] - centroid;
            if (diff.x*diff.x+diff.y*diff.y+diff.z*diff.z < circum_r*circum_r*4.0)
                mean_H_local += orig_H[v];
        }
        cf.mean_H_local = mean_H_local;

        // ---------------------------------------------------------------
        // Pottmann Eq.7 — solve for sphere centre C.
        // The 7 sample points and their normals constrain C as the point
        // where all normals are orthogonal to (v - C).
        //
        // After solving, sign is fixed using the mean outward normal:
        // C must be on the INWARD side → dot(C - centroid, mean_out) < 0.
        // ---------------------------------------------------------------

        // Compute mean outward normal from the 7 samples (reliably outward)
        GEO::vec3 mean_out(0,0,0);
        int n = (int)nbr_pts.size();
        for (int i = 0; i < n; ++i) mean_out = mean_out + nbr_nrm[i];
        double ml = GEO::length(mean_out);
        if (ml > 1e-14) mean_out = mean_out * (1.0/ml);
        else mean_out = coarse_normal;

        // mean H for descriptor
        double mean_H = 0.0; int H_count = 0;
        for (GEO::index_t v = 0; v < nv_orig; ++v) {
            GEO::vec3 diff = orig_pts[v] - centroid;
            if (diff.x*diff.x+diff.y*diff.y+diff.z*diff.z < circum_r*circum_r*9.0)
                { mean_H += orig_H[v]; ++H_count; }
        }
        if (H_count > 0) mean_H /= H_count;
        cf.mean_H_local = mean_H;

        // Build and solve 3×3 linear system (Pottmann Eq.7)
        double A[3][3] = {};
        double b3[3]   = {};
        for (int i = 0; i < n; ++i) {
            const GEO::vec3& ni = nbr_nrm[i];
            const GEO::vec3& vi = nbr_pts[i];
            double nn[3] = {ni.x, ni.y, ni.z};
            for (int r = 0; r < 3; ++r)
                for (int c2 = 0; c2 < 3; ++c2)
                    A[r][c2] += (r==c2 ? 1.0 : 0.0) - nn[r]*nn[c2];
            double ndotv = ni.x*vi.x + ni.y*vi.y + ni.z*vi.z;
            for (int r = 0; r < 3; ++r)
                b3[r] += (r==0?vi.x:r==1?vi.y:vi.z) - nn[r]*ndotv;
        }
        double Msys3[3][4];
        for (int r = 0; r < 3; ++r) {
            for (int c2 = 0; c2 < 3; ++c2) Msys3[r][c2] = A[r][c2];
            Msys3[r][3] = b3[r];
        }
        bool ok3 = true;
        for (int col = 0; col < 3 && ok3; ++col) {
            int pivot = col;
            for (int row = col+1; row < 3; ++row)
                if (std::fabs(Msys3[row][col]) > std::fabs(Msys3[pivot][col]))
                    pivot = row;
            if (std::fabs(Msys3[pivot][col]) < 1e-12) { ok3 = false; break; }
            std::swap(Msys3[col], Msys3[pivot]);
            double inv = 1.0 / Msys3[col][col];
            for (int j = col; j <= 3; ++j) Msys3[col][j] *= inv;
            for (int row = 0; row < 3; ++row) {
                if (row == col) continue;
                double fac = Msys3[row][col];
                for (int j = col; j <= 3; ++j) Msys3[row][j] -= fac*Msys3[col][j];
            }
        }
        if (!ok3) continue;
        GEO::vec3 C(Msys3[0][3], Msys3[1][3], Msys3[2][3]);

        // No sign flip needed: Pottmann Eq.7 places C correctly when normals
        // are consistently oriented. The tessellation projects pts away from C
        // regardless of which side C is on (d = pt - C, then C + d*(R/len)).
        // If C ends up on the wrong side, the face will appear inverted —
        // which is a mesh winding issue, not a sphere fitting issue.

        // R = mean distance from sample points to C
        double sum_r = 0.0;
        for (auto& q : nbr_pts) sum_r += GEO::length(q - C);
        double R = sum_r / (double)n;

        // Reject if R is implausibly small
        if (R < 1e-10) continue;

        cf.sphere_center   = C;
        cf.sphere_radius   = R;
        cf.valid           = true;
        cf.sphere_residual = 0.0;
    }
}

// =======================================================================
// JOINT SPHERE OPTIMISATION (Pottmann 2023, Sec 4)
// =======================================================================
//
// After independent per-face sphere fitting (Pottmann Eq.7), adjacent
// spheres are independent and do NOT necessarily intersect in a real
// circle → shared arc construction is ill-defined → tessellation has
// cracks/folds.
//
// Fix: joint optimisation that moves each sphere's centre C_f and
// radius R_f so that:
//   E_sphere: coarse vertices lie on their face's sphere
//   E_smooth: adjacent spheres intersect at a smooth angle (angle → 0)
//
// We use alternating least-squares (Gauss-Newton style):
//   Step A: fix coarse vertex positions, update sphere (C,R) per face
//           by least-squares through the 3 vertices.
//   Step B: fix spheres, move each coarse vertex to its mean projection
//           across all faces that share it.
// Repeat for max_iter iterations. This converges when coarse vertices
// lie on all their adjacent spheres simultaneously.
//
// After convergence, every shared edge is the intersection circle of
// the two adjacent spheres — the shared arc is well-defined and lies
// on BOTH spheres exactly.

void optimise_joint_spheres(
    GEO::Mesh&                 coarse_mesh,
    std::vector<CurvedFace>&   cfaces,
    int                        max_iter)
{
    const GEO::index_t nf = coarse_mesh.facets.nb();
    const GEO::index_t nv = coarse_mesh.vertices.nb();

    // Working copy of coarse vertex positions
    std::vector<GEO::vec3> vpos(nv);
    for (GEO::index_t v = 0; v < nv; ++v)
        vpos[v] = vertex_point(coarse_mesh, v);

    // Build vertex → faces map (only valid/interior faces)
    std::vector<std::vector<GEO::index_t>> v2f(nv);
    for (GEO::index_t f = 0; f < nf; ++f) {
        if (!cfaces[f].valid) continue;
        for (int k = 0; k < 3; ++k)
            v2f[cfaces[f].v[k]].push_back(f);
    }

    // Mark boundary vertices — on any edge with only 1 adjacent face.
    // These vertices must NOT move in Step B: they anchor the open boundary.
    std::vector<bool> is_boundary_v(nv, false);
    for (GEO::index_t f = 0; f < nf; ++f) {
        if (coarse_mesh.facets.nb_vertices(f) < 3) continue;
        for (GEO::index_t le = 0; le < coarse_mesh.facets.nb_vertices(f); ++le) {
            if (coarse_mesh.facets.adjacent(f, le) == GEO::NO_FACET) {
                // This edge is a boundary edge — mark both its vertices
                is_boundary_v[coarse_mesh.facets.vertex(f, le)] = true;
                is_boundary_v[coarse_mesh.facets.vertex(f, (le+1)%3)] = true;
            }
        }
    }

    for (int iter = 0; iter < max_iter; ++iter) {

        // ---- Step A: refit sphere to current 3 coarse vertices ----
        // For a triangle with 3 points + their original normals, solve
        // the minimal-residual sphere using 4×4 LS (sphere through 3
        // pts is underdetermined in R³ — we add the centroid normal
        // constraint to pin the solution). Simpler: use the closed-form
        // circumsphere lifted by mean curvature, re-evaluated at current
        // vertex positions.
        for (GEO::index_t f = 0; f < nf; ++f) {
            CurvedFace& cf = cfaces[f];
            if (!cf.valid) continue;

            GEO::vec3 p[3];
            for (int k = 0; k < 3; ++k) p[k] = vpos[cf.v[k]];

            // Fit sphere through exactly 3 points + normal direction.
            // The 3 points define a circle; the sphere centre lies on the
            // line through the circumcenter perpendicular to the triangle.
            // We choose the sphere centre as:
            //   C = circumcenter + t * face_normal
            // where t minimises |C - p[i]|² consistency.
            // For 3 points on a sphere: |C-p0|=|C-p1|=|C-p2|=R.
            // The circumcenter already satisfies |cc-p0|=|cc-p1|=|cc-p2|=rcc.
            // With C = cc + t*n: R² = rcc² + t².
            // t is a free parameter — we choose t to minimise displacement
            // from the previous sphere centre (stability).
            GEO::vec3 e1 = p[1]-p[0], e2 = p[2]-p[0];
            GEO::vec3 n  = GEO::cross(e1, e2);
            double n_len = GEO::length(n);
            if (n_len < 1e-14) continue;
            GEO::vec3 n_unit = n * (1.0/n_len);

            // Circumcenter in the triangle plane
            double d1 = GEO::dot(e1,e1), d2 = GEO::dot(e2,e2), d12 = GEO::dot(e1,e2);
            double denom = 2.0 * n_len * n_len;
            double l1 = d2*(d1-d12)/denom, l2 = d1*(d2-d12)/denom;
            GEO::vec3 cc = p[0] + e1*l1 + e2*l2;
            double rcc = GEO::length(cc - p[0]);

            // Choose t so C is closest to previous sphere centre
            GEO::vec3 prev_C = cf.sphere_center;
            double t = GEO::dot(prev_C - cc, n_unit);

            // Clamp t to keep the sphere from becoming degenerate
            // Min |t| ensures sphere has visible curvature
            // (don't let it collapse to circumscribed circle)
            double t_min_abs = rcc * 0.05; // at least 5% lift
            if (std::fabs(t) < t_min_abs)
                t = (t >= 0.0) ? t_min_abs : -t_min_abs;

            GEO::vec3 C = cc + n_unit * t;
            double R = std::sqrt(rcc*rcc + t*t);

            cf.sphere_center = C;
            cf.sphere_radius  = R;
            cf.sphere_residual = 0.0;
        }

        // ---- Step B: move each coarse vertex to mean projection ----
        // Each vertex is shared by several faces, each with a different sphere.
        // Move the vertex to the mean of its projections onto all adjacent
        // spheres. This pulls the vertex toward the intersection of all adjacent
        // spheres — the converged position lies on all of them simultaneously.
        std::vector<GEO::vec3> new_vpos = vpos;
        for (GEO::index_t v = 0; v < nv; ++v) {
            if (v2f[v].empty()) continue;
            if (is_boundary_v[v]) continue;  // never move boundary vertices
            GEO::vec3 mean_proj(0,0,0);
            int cnt = 0;
            for (GEO::index_t f : v2f[v]) {
                const CurvedFace& cf = cfaces[f];
                if (!cf.valid) continue;
                GEO::vec3 d = vpos[v] - cf.sphere_center;
                double len = GEO::length(d);
                if (len > 1e-14) {
                    mean_proj = mean_proj + cf.sphere_center + d*(cf.sphere_radius/len);
                    ++cnt;
                }
            }
            if (cnt > 0)
                new_vpos[v] = mean_proj * (1.0/cnt);
        }

        // Compute max vertex displacement for convergence check
        double max_disp = 0.0;
        for (GEO::index_t v = 0; v < nv; ++v)
            max_disp = std::max(max_disp, GEO::length(new_vpos[v] - vpos[v]));

        vpos = new_vpos;

        if (max_disp < 1e-8) {
            GEO::Logger::out("CurveKID")
                << "  Joint opt converged at iter " << iter
                << " (max_disp=" << max_disp << ")" << std::endl;
            break;
        }
    }

    // Write optimised positions back to coarse_mesh
    for (GEO::index_t v = 0; v < nv; ++v) {
        double* p = coarse_mesh.vertices.point_ptr(v);
        p[0] = vpos[v].x; p[1] = vpos[v].y; p[2] = vpos[v].z;
    }

    // Update cfaces.v positions (they reference the coarse mesh)
    // The face vertex indices stay the same — positions updated above.

    GEO::Logger::out("CurveKID")
        << "  Joint sphere optimisation done ("
        << max_iter << " max iters)" << std::endl;
}

// =======================================================================
// CURVED MESH CONSTRUCTION — Stage 3: tessellate onto sphere for viz
// =======================================================================

void tessellate_curved_faces(
    const GEO::Mesh&               coarse_mesh,
    const std::vector<CurvedFace>& cfaces,
    GEO::Mesh&                     viz_mesh,
    int                            subdiv)
{
    viz_mesh.clear();
    if (subdiv < 1) subdiv = 1;

    // Threshold: if sphere radius > this multiple of the coarse triangle's
    // circumradius the face is effectively flat — skip sphere projection
    // and use the flat triangle directly (avoids degenerate huge-sphere output).
    const double FLAT_THRESHOLD = 500.0;

    GEO::Attribute<int> attr_face_id(viz_mesh.facets.attributes(), "curved_face_id");

    for (GEO::index_t f = 0; f < (GEO::index_t)cfaces.size(); ++f) {
        const CurvedFace& cf = cfaces[f];
        if (!cf.valid) continue;
        if (coarse_mesh.facets.nb_vertices(cf.face_id) < 3) continue;

        // Coarse triangle corners
        GEO::vec3 p[3];
        for (int k = 0; k < 3; ++k)
            p[k] = vertex_point(coarse_mesh, cf.v[k]);

        GEO::vec3 C = cf.sphere_center;
        double    R = cf.sphere_radius;

        // Circumradius of coarse triangle (used to detect degenerate spheres)
        double ea = GEO::length(p[1]-p[0]);
        double eb = GEO::length(p[2]-p[1]);
        double ec = GEO::length(p[0]-p[2]);
        double s  = (ea+eb+ec)*0.5;
        double area_tri = std::sqrt(std::max(0.0, s*(s-ea)*(s-eb)*(s-ec)));
        double circum_r = (area_tri > 1e-14) ? (ea*eb*ec)/(4.0*area_tri)
                                              : (ea+eb+ec)/3.0;

        bool use_flat = (R > FLAT_THRESHOLD * circum_r) || (R < 1e-10);

        // Project a point onto the fitted sphere (radial projection from C)
        // For flat faces, identity projection.
        auto proj = [&](GEO::vec3 pt) -> GEO::vec3 {
            if (use_flat) return pt;
            GEO::vec3 d = pt - C;
            double len = GEO::length(d);
            if (len < 1e-14) return pt;
            return C + d * (R / len);
        };

        // Build barycentric grid: interpolate in FLAT triangle space, then project.
        // This is correct for any sphere radius and avoids blending on-sphere
        // positions (which degenerates when C is far away).
        int n_verts = (subdiv+1)*(subdiv+2)/2;
        GEO::index_t base = viz_mesh.vertices.create_vertices(n_verts);

        int idx = 0;
        for (int i = 0; i <= subdiv; ++i) {
            for (int j = 0; j <= i; ++j) {
                double bu = (subdiv > 0) ? double(i-j)/double(subdiv) : 0.0;
                double bv = (subdiv > 0) ? double(j)  /double(subdiv) : 0.0;
                double bw = 1.0 - bu - bv;
                // Interpolate in flat triangle, then radially project onto sphere
                GEO::vec3 pt = p[0]*bw + p[1]*bu + p[2]*bv;
                pt = proj(pt);
                double* vp = viz_mesh.vertices.point_ptr(base + idx);
                vp[0] = pt.x; vp[1] = pt.y; vp[2] = pt.z;
                ++idx;
            }
        }

        // Grid index: row i, col j → i*(i+1)/2 + j
        // Triangulate
        for (int i = 0; i < subdiv; ++i) {
            for (int j = 0; j <= i; ++j) {
                int g00 = i*(i+1)/2 + j;
                int g10 = (i+1)*(i+2)/2 + j;
                int g11 = (i+1)*(i+2)/2 + j+1;
                // Upward triangle
                GEO::index_t nf = viz_mesh.facets.create_triangle(
                    base + g00, base + g10, base + g11);
                attr_face_id[nf] = (int)f;
                // Downward triangle
                if (j < i) {
                    int g01 = i*(i+1)/2 + j+1;
                    GEO::index_t nf2 = viz_mesh.facets.create_triangle(
                        base + g00, base + g11, base + g01);
                    attr_face_id[nf2] = (int)f;
                }
            }
        }
    }

    viz_mesh.facets.connect();
}

// =======================================================================
// CURVED MESH CONSTRUCTION — Stage 4: approximation error
// =======================================================================

void compute_approximation_error(
    const GEO::Mesh&               original,
    const std::vector<CurvedFace>& cfaces,
    std::vector<double>&           errors)
{
    const GEO::index_t nv = original.vertices.nb();
    errors.assign(nv, 0.0);

    // For each original vertex, find nearest curved face sphere
    // and compute distance to that sphere surface.
    // Simple O(nv * nf) scan — fast enough for typical mesh sizes.
    for (GEO::index_t v = 0; v < nv; ++v) {
        GEO::vec3 pt = vertex_point(original, v);
        double best_err = std::numeric_limits<double>::max();
        for (auto& cf : cfaces) {
            if (!cf.valid) continue;
            double dist_to_sphere =
                std::fabs(GEO::length(pt - cf.sphere_center) - cf.sphere_radius);
            best_err = std::min(best_err, dist_to_sphere);
        }
        errors[v] = (best_err < std::numeric_limits<double>::max()) ? best_err : 0.0;
    }
}

// =======================================================================
// PLANAR APPROXIMATION ERROR
// For each original vertex, compute distance to the nearest triangle
// of the coarse flat mesh.  This is E_planar(N) in the paper.
// =======================================================================

static double point_to_triangle_dist(
    const GEO::vec3& p,
    const GEO::vec3& a,
    const GEO::vec3& b,
    const GEO::vec3& c)
{
    // Project p onto the plane of (a,b,c), clamp to triangle, return distance.
    GEO::vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = GEO::dot(ab, ap), d2 = GEO::dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return GEO::length(p - a);

    GEO::vec3 bp = p - b;
    double d3 = GEO::dot(ab, bp), d4 = GEO::dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return GEO::length(p - b);

    GEO::vec3 cp = p - c;
    double d5 = GEO::dot(ab, cp), d6 = GEO::dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return GEO::length(p - c);

    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        return GEO::length(p - (a + ab * v));
    }
    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        return GEO::length(p - (a + ac * w));
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
        double w = (d4-d3) / ((d4-d3) + (d5-d6));
        return GEO::length(p - (b + (c - b) * w));
    }
    // Inside triangle: distance to plane
    double denom = va + vb + vc;
    if (std::fabs(denom) < 1e-14) return 0.0;
    double ww = vc / denom;
    double wv = vb / denom;
    GEO::vec3 closest = a + ab * ww + ac * wv;
    return GEO::length(p - closest);
}

void compute_planar_approximation_error(
    const GEO::Mesh& original,
    const GEO::Mesh& coarse_flat,
    std::vector<double>& errors)
{
    const GEO::index_t nv = original.vertices.nb();
    const GEO::index_t nf = coarse_flat.facets.nb();
    errors.assign(nv, 0.0);

    for (GEO::index_t v = 0; v < nv; ++v) {
        GEO::vec3 pt = vertex_point(original, v);
        double best = std::numeric_limits<double>::max();
        for (GEO::index_t f = 0; f < nf; ++f) {
            if (coarse_flat.facets.nb_vertices(f) < 3) continue;
            GEO::vec3 a = vertex_point(coarse_flat, coarse_flat.facets.vertex(f,0));
            GEO::vec3 b = vertex_point(coarse_flat, coarse_flat.facets.vertex(f,1));
            GEO::vec3 c = vertex_point(coarse_flat, coarse_flat.facets.vertex(f,2));
            best = std::min(best, point_to_triangle_dist(pt, a, b, c));
        }
        errors[v] = (best < std::numeric_limits<double>::max()) ? best : 0.0;
    }
}

// =======================================================================
// CIRCUMSPHERE through 3 vertices with normal offset
// =======================================================================

bool fit_circumsphere(
    const GEO::vec3& p0,
    const GEO::vec3& p1,
    const GEO::vec3& p2,
    double           t,
    GEO::vec3&       center,
    double&          radius)
{
    // Triangle edge vectors
    GEO::vec3 e1 = p1 - p0;
    GEO::vec3 e2 = p2 - p0;
    GEO::vec3 n  = GEO::cross(e1, e2);
    double n_len = GEO::length(n);
    if (n_len < 1e-14) return false;   // degenerate triangle

    GEO::vec3 n_unit = n * (1.0 / n_len);

    // Circumcenter of the triangle (in the plane of the triangle).
    // Formula: circumcenter = p0 + (|e2|^2*(e1·e1-e1·e2)*e2 - |e1|^2*(e2·e2-e1·e2)*e1) / (2*|n|^2)
    // Simplified using cross products:
    double d1 = GEO::dot(e1, e1);
    double d2 = GEO::dot(e2, e2);
    double d12= GEO::dot(e1, e2);
    // circumcenter via barycentric approach
    // lambda1 = d2*(d1-d12) / (2*n_len^2)
    // lambda2 = d1*(d2-d12) / (2*n_len^2)
    double denom = 2.0 * n_len * n_len;
    double l1 = d2 * (d1 - d12) / denom;
    double l2 = d1 * (d2 - d12) / denom;
    GEO::vec3 cc = p0 + e1 * l1 + e2 * l2;  // circumcenter in plane

    // Circumradius in plane
    double cc_r = GEO::length(cc - p0);

    // Sphere center = circumcenter + t * circumradius * face_normal
    // t controls the "lift" off the plane.
    // R^2 = circumradius^2 + (t*circumradius)^2 = circumradius^2*(1+t^2)
    double lift = t * cc_r;
    center = cc + n_unit * lift;
    radius = std::sqrt(cc_r * cc_r + lift * lift);

    // Verify: all 3 vertices should be on the sphere
    // (they are by construction — circumcenter equidistant to p0,p1,p2)
    return (radius > 1e-14);
}

// =======================================================================
// SYNTHETIC MODE: fit circumsphere per coarse face, preserve exact corners
// =======================================================================

void fit_curved_faces_synthetic(
    const GEO::Mesh&         coarse_mesh,
    std::vector<CurvedFace>& cfaces,
    double                   t_factor)
{
    const GEO::index_t nf = coarse_mesh.facets.nb();
    cfaces.resize(nf);

    for (GEO::index_t f = 0; f < nf; ++f) {
        CurvedFace& cf = cfaces[f];
        cf.face_id = f;
        cf.valid   = false;

        if (coarse_mesh.facets.nb_vertices(f) < 3) continue;

        cf.v[0] = coarse_mesh.facets.vertex(f, 0);
        cf.v[1] = coarse_mesh.facets.vertex(f, 1);
        cf.v[2] = coarse_mesh.facets.vertex(f, 2);

        GEO::vec3 p0 = vertex_point(coarse_mesh, cf.v[0]);
        GEO::vec3 p1 = vertex_point(coarse_mesh, cf.v[1]);
        GEO::vec3 p2 = vertex_point(coarse_mesh, cf.v[2]);

        GEO::vec3 C; double R;
        if (fit_circumsphere(p0, p1, p2, t_factor, C, R)) {
            cf.sphere_center   = C;
            cf.sphere_radius   = R;
            // Residual = max deviation of corners from sphere (should be ~0)
            double e0 = std::fabs(GEO::length(p0-C) - R);
            double e1 = std::fabs(GEO::length(p1-C) - R);
            double e2 = std::fabs(GEO::length(p2-C) - R);
            cf.sphere_residual = std::max({e0, e1, e2});
            cf.valid = true;
        }
    }
}

// =======================================================================
// SHARED ARC RESOLUTION
// =======================================================================

void resolve_shared_arcs(
    const GEO::Mesh&               coarse_mesh,
    const std::vector<CurvedFace>& cfaces,
    std::vector<GEO::index_t>&     edge_winner)
{
    const GEO::index_t nf = coarse_mesh.facets.nb();
    // Initialise: each edge's winner is the face itself
    edge_winner.assign((size_t)nf * 3u, GEO::index_t(-1));
    for (GEO::index_t f = 0; f < nf; ++f)
        for (GEO::index_t e = 0; e < 3; ++e)
            edge_winner[(size_t)f*3u+e] = f;

    // Build edge → (face, local_edge) map using sorted vertex pair as key
    struct PairHash {
        size_t operator()(std::pair<GEO::index_t,GEO::index_t> p) const {
            return std::hash<uint64_t>()(
                (uint64_t)p.first << 32 | (uint64_t)p.second);
        }
    };
    // Map: canonical edge → (face0, local_edge0)
    std::unordered_map<
        std::pair<GEO::index_t,GEO::index_t>,
        std::pair<GEO::index_t,GEO::index_t>,
        PairHash> edge_map;

    for (GEO::index_t f = 0; f < nf; ++f) {
        if (coarse_mesh.facets.nb_vertices(f) < 3) continue;
        for (GEO::index_t e = 0; e < 3; ++e) {
            GEO::index_t va = coarse_mesh.facets.vertex(f, e);
            GEO::index_t vb = coarse_mesh.facets.vertex(f, (e+1)%3);
            if (va > vb) std::swap(va, vb);
            auto key = std::make_pair(va, vb);
            auto it  = edge_map.find(key);
            if (it == edge_map.end()) {
                edge_map[key] = std::make_pair(f, e);
            } else {
                // Second face sharing this edge: compare residuals
                GEO::index_t f0 = it->second.first;
                GEO::index_t e0 = it->second.second;
                // Winner = lower sphere_residual (better fit)
                GEO::index_t winner = f;
                if (cfaces[f0].valid && cfaces[f].valid) {
                    winner = (cfaces[f0].sphere_residual <= cfaces[f].sphere_residual)
                             ? f0 : f;
                } else if (cfaces[f0].valid) {
                    winner = f0;
                }
                edge_winner[(size_t)f *3u + e]  = winner;
                edge_winner[(size_t)f0*3u + e0] = winner;
            }
        }
    }
}

// =======================================================================
// TESSELLATION v2: shared-arc aware
// =======================================================================

void tessellate_curved_faces_v2(
    const GEO::Mesh&                   coarse_mesh,
    const std::vector<CurvedFace>&     cfaces,
    const std::vector<GEO::index_t>&   /*edge_winner — unused, shared-edge approach*/,
    GEO::Mesh&                         viz_mesh,
    int                                subdiv)
{
    // -----------------------------------------------------------------------
    // WATERTIGHT shared-edge tessellation (same approach as build_sphere_ground_truth).
    //
    // Key insight: cracks come from each face creating its OWN edge vertices
    // independently. Fix: generate edge interior points ONCE and share them
    // between adjacent faces. This guarantees identical vertex positions on
    // shared edges → zero cracks, zero holes.
    //
    // For each coarse edge (va,vb):
    //   - Use the LOWER-index face's sphere (or the face that has a valid sphere)
    //   - Generate (subdiv-1) slerp points on that sphere
    //   - Both adjacent faces look up the same points
    //
    // Corners are shared trivially — they map to the same coarse vertex index.
    // Interior points are face-private (each face uses its own sphere there).
    // -----------------------------------------------------------------------
    viz_mesh.clear();
    if (subdiv < 1) subdiv = 1;

    GEO::Attribute<int> attr_face_id(viz_mesh.facets.attributes(), "curved_face_id");
    GEO::Attribute<int> attr_is_curved(viz_mesh.facets.attributes(), "is_curved");

    const double FLAT_THRESHOLD = 5000.0;

    // Precompute circumradius per coarse face (for flat-fallback threshold)
    std::vector<double> coarse_circum_r(cfaces.size(), 0.0);
    for (GEO::index_t f = 0; f < (GEO::index_t)cfaces.size(); ++f) {
        const CurvedFace& cf = cfaces[f];
        GEO::vec3 p0 = vertex_point(coarse_mesh, cf.v[0]);
        GEO::vec3 p1 = vertex_point(coarse_mesh, cf.v[1]);
        GEO::vec3 p2 = vertex_point(coarse_mesh, cf.v[2]);
        double a = GEO::length(p1-p0), b = GEO::length(p2-p1), c2 = GEO::length(p0-p2);
        double s  = (a+b+c2)*0.5;
        double ar = std::sqrt(std::max(0.0, s*(s-a)*(s-b)*(s-c2)));
        coarse_circum_r[f] = (ar > 1e-14) ? (a*b*c2)/(4.0*ar) : (a+b+c2)/3.0;
    }

    // Helper: project point onto a face's sphere (with flat fallback)
    auto proj_onto_face = [&](GEO::index_t fw, GEO::vec3 pt) -> GEO::vec3 {
        if (fw >= (GEO::index_t)cfaces.size() || !cfaces[fw].valid) return pt;
        const CurvedFace& cf = cfaces[fw];
        if (cf.sphere_radius > FLAT_THRESHOLD * coarse_circum_r[fw]) return pt;
        GEO::vec3 d = pt - cf.sphere_center;
        double len = GEO::length(d);
        if (len < 1e-14) return pt;
        return cf.sphere_center + d * (cf.sphere_radius / len);
    };

    // Helper: slerp on a face's sphere between two points
    auto slerp_on_face = [&](GEO::index_t fw, GEO::vec3 pa, GEO::vec3 pb, double t) -> GEO::vec3 {
        if (fw >= (GEO::index_t)cfaces.size() || !cfaces[fw].valid)
            return pa*(1-t) + pb*t;
        const CurvedFace& cf = cfaces[fw];
        // If sphere is flat (very large R relative to face size), use linear interp
        if (coarse_circum_r[fw] > 1e-14 && cf.sphere_radius > FLAT_THRESHOLD * coarse_circum_r[fw])
            return pa*(1-t) + pb*t;
        GEO::vec3 C = cf.sphere_center;
        double R    = cf.sphere_radius;
        GEO::vec3 da = pa - C, db = pb - C;
        double la = GEO::length(da), lb = GEO::length(db);
        if (la < 1e-14 || lb < 1e-14) return pa*(1-t)+pb*t;
        GEO::vec3 ua = da*(1.0/la), ub = db*(1.0/lb);
        double dot = std::max(-1.0, std::min(1.0, GEO::dot(ua, ub)));
        double theta = std::acos(dot);
        if (std::fabs(theta) < 1e-10) return C + ua*R;
        GEO::vec3 interp = ua*(std::sin((1-t)*theta)/std::sin(theta))
                         + ub*(std::sin(t*theta)    /std::sin(theta));
        return C + interp * R;
    };

    // --- Step 1: add coarse vertices to viz_mesh (index == coarse vertex index) ---
    // Corner vertices are shared by multiple faces. We store them at their flat
    // coarse position first; they will be projected face-by-face in Step 3.
    // BUT since multiple faces share a corner, we cannot project each corner
    // onto a different sphere — we keep corners at their flat coarse positions.
    // This is acceptable: the coarse vertices already lie close to the original
    // surface (they come from remeshing), so the error at corners is small.
    const GEO::index_t nv_c = coarse_mesh.vertices.nb();
    viz_mesh.vertices.create_vertices(nv_c);
    for (GEO::index_t v = 0; v < nv_c; ++v) {
        const double* p = coarse_mesh.vertices.point_ptr(v);
        double* vp = viz_mesh.vertices.point_ptr(v);
        vp[0]=p[0]; vp[1]=p[1]; vp[2]=p[2];
    }

    // --- Step 2: build edge point table ---
    // For coarse edge (va,vb) with va<vb: (subdiv-1) interior points
    // stored in direction va→vb, projected onto the sphere of the
    // lower-index valid curved face that owns this edge.
    //
    // "owns" = the first face (lower face index) among all cfaces that has
    // this edge AND has a valid sphere.  If neither is valid, use flat lerp.

    // Build edge→face map: for each directed edge (u,v), which cface index?
    std::map<std::pair<GEO::index_t,GEO::index_t>, GEO::index_t> edge_to_face;
    for (GEO::index_t f = 0; f < (GEO::index_t)cfaces.size(); ++f) {
        const CurvedFace& cf = cfaces[f];
        for (int e = 0; e < 3; ++e) {
            GEO::index_t va = cf.v[e], vb = cf.v[(e+1)%3];
            // Store canonical (min,max) → owner face (prefer lower index with valid sphere)
            GEO::index_t lo = std::min(va,vb), hi = std::max(va,vb);
            auto key = std::make_pair(lo, hi);
            auto it = edge_to_face.find(key);
            if (it == edge_to_face.end()) {
                edge_to_face[key] = f;
            } else {
                // Keep whichever face has a valid sphere; on tie keep lower index
                if (!cfaces[it->second].valid && cfaces[f].valid)
                    it->second = f;
            }
        }
    }

    // edge_pts: canonical edge → vector of (subdiv-1) viz_mesh vertex indices, dir lo→hi
    // ALL edge points are projected onto the OWNER face's sphere — same as interior.
    // This means ALL points on a face patch lie on its sphere → no tent distortion.
    // Adjacent faces that share an edge use THE SAME edge points (same owner sphere)
    // → watertight by construction.
    std::map<std::pair<GEO::index_t,GEO::index_t>, std::vector<GEO::index_t>> edge_pts;

    auto get_edge_pts = [&](GEO::index_t va, GEO::index_t vb)
        -> std::vector<GEO::index_t>
    {
        bool flipped = (va > vb);
        GEO::index_t lo = flipped ? vb : va;
        GEO::index_t hi = flipped ? va : vb;
        auto key = std::make_pair(lo, hi);
        auto it = edge_pts.find(key);
        if (it == edge_pts.end()) {
            // Find owner face — the first (lowest-index) valid face owning this edge
            GEO::index_t owner = GEO::NO_FACET;
            auto fit = edge_to_face.find(key);
            if (fit != edge_to_face.end()) owner = fit->second;

            GEO::vec3 plo = vertex_point(coarse_mesh, lo);
            GEO::vec3 phi_v = vertex_point(coarse_mesh, hi);

            std::vector<GEO::index_t> ids;
            for (int s = 1; s < subdiv; ++s) {
                double t = double(s) / double(subdiv);
                GEO::vec3 pt;
                if (owner != GEO::NO_FACET && cfaces[owner].valid) {
                    // After joint optimisation, plo and phi_v both lie on the
                    // owner sphere. Use slerp on the owner sphere → arc lies
                    // exactly on the sphere surface. This is the correct shared
                    // circular arc from Pottmann 2023.
                    pt = slerp_on_face(owner, plo, phi_v, t);
                } else {
                    pt = plo*(1.0-t) + phi_v*t;
                }
                GEO::index_t idx = viz_mesh.vertices.create_vertices(1);
                double* vp = viz_mesh.vertices.point_ptr(idx);
                vp[0]=pt.x; vp[1]=pt.y; vp[2]=pt.z;
                ids.push_back(idx);
            }
            edge_pts[key] = ids;
            it = edge_pts.find(key);
        }
        std::vector<GEO::index_t> result = it->second;
        if (flipped) std::reverse(result.begin(), result.end());
        return result;
    };

    // --- Step 3: tessellate each face using shared edges ---
    int S = subdiv;
    for (GEO::index_t f = 0; f < (GEO::index_t)cfaces.size(); ++f) {
        const CurvedFace& cf = cfaces[f];
        if (coarse_mesh.facets.nb_vertices(cf.face_id) < 3) continue;

        // Grid: grid[i][j] where i+j <= S
        // i=0,j=0 → cv[0];  i=0,j=S → cv[1];  i=S,j=0 → cv[2]
        std::vector<std::vector<GEO::index_t>> grid(S+1,
            std::vector<GEO::index_t>(S+1, GEO::NO_VERTEX));

        grid[0][0] = cf.v[0];
        grid[0][S] = cf.v[1];
        grid[S][0] = cf.v[2];

        // Edge cv[0]→cv[1] (i=0, j=1..S-1)
        { auto ep = get_edge_pts(cf.v[0], cf.v[1]);
          for (int j=1;j<S;++j) grid[0][j] = ep[j-1]; }

        // Edge cv[0]→cv[2] (j=0, i=1..S-1)
        { auto ep = get_edge_pts(cf.v[0], cf.v[2]);
          for (int i=1;i<S;++i) grid[i][0] = ep[i-1]; }

        // Edge cv[1]→cv[2] (i+j=S, i=1..S-1)
        { auto ep = get_edge_pts(cf.v[1], cf.v[2]);
          for (int i=1;i<S;++i) grid[i][S-i] = ep[i-1]; }

        // Interior points on this face's own sphere
        for (int i=1;i<S;++i) {
            for (int j=1;j<S-i;++j) {
                double bv0 = 1.0 - double(i+j)/S;
                double bv1 = double(j)/S;
                double bv2 = double(i)/S;
                GEO::vec3 p0 = vertex_point(coarse_mesh, cf.v[0]);
                GEO::vec3 p1 = vertex_point(coarse_mesh, cf.v[1]);
                GEO::vec3 p2 = vertex_point(coarse_mesh, cf.v[2]);
                GEO::vec3 pt_flat = p0*bv0 + p1*bv1 + p2*bv2;
                GEO::vec3 pt = proj_onto_face(f, pt_flat);
                GEO::index_t idx = viz_mesh.vertices.create_vertices(1);
                double* vp = viz_mesh.vertices.point_ptr(idx);
                vp[0]=pt.x; vp[1]=pt.y; vp[2]=pt.z;
                grid[i][j] = idx;
            }
        }

        // Build sub-triangles
        int curved_flag = cf.valid ? 1 : 0;
        for (int i=0; i<S; ++i) {
            for (int j=0; j<S-i; ++j) {
                GEO::index_t a=grid[i][j], b=grid[i][j+1], c2=grid[i+1][j];
                if (a!=GEO::NO_VERTEX && b!=GEO::NO_VERTEX && c2!=GEO::NO_VERTEX) {
                    GEO::index_t nf2 = viz_mesh.facets.create_triangle(a, b, c2);
                    attr_face_id[nf2]   = (int)f;
                    attr_is_curved[nf2] = curved_flag;
                }
                if (i+j+2 <= S) {
                    GEO::index_t d2=grid[i+1][j+1];
                    if (b!=GEO::NO_VERTEX && d2!=GEO::NO_VERTEX && c2!=GEO::NO_VERTEX) {
                        GEO::index_t nf3 = viz_mesh.facets.create_triangle(b, d2, c2);
                        attr_face_id[nf3]   = (int)f;
                        attr_is_curved[nf3] = curved_flag;
                    }
                }
            }
        }
    }

    viz_mesh.facets.connect();
}

} // namespace CurveKIDAlgo
} // namespace OGF
