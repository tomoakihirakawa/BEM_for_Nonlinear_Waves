#ifndef GNUPLOT_H
   #define GNUPLOT_H

   #include <stdio.h>
   #include <stdlib.h>
   #include <iostream>
   #include <map>
   #include <sstream>
   #include <string>
   #include <vector>

   #ifndef _USE_MATH_DEFINES
      #define _USE_MATH_DEFINES
      #include <cmath>
      #ifndef M_PI
         #define M_PI (3.14159265358979323846)
      #endif
   #endif

class GNUPLOT {
  public:
   FILE* proc;
   int counter, nplot, n_decld_plot;
   std::string decld_str;
   std::vector<std::vector<std::vector<double>>> Data;  // (3){(2){(1)vecx,vecy,vecz}}
   std::vector<std::map<std::string, std::string>> OPT;
   std::vector<std::vector<std::vector<std::vector<double>>>> SplotData;
   // ex.
   // PLOT.SaveSplotData( { {1, 1, 2}, {2, 1, 2}, {3, 1, 2}, {4, 1, 1} },
   //                       {1, 2, 2}, {2, 2, 2}, {3, 2, 2}, {4, 2, 2} },
   //                       {1, 3, 2}, {2, 3, 2}, {3, 3, 2}, {4, 3, 1} } })
   std::vector<std::map<std::string, std::string>> SplotOPT;
   std::vector<std::vector<std::vector<std::vector<double>>>> VectorData;  // (4){(3){(2){(1)vecx1, vecy1},{vecx2, vecy2}}}
   std::vector<std::map<std::string, std::string>> VectorOPT;
   std::vector<std::vector<std::vector<double>>> MatrixData;
   // ex.
   // PLOT.SaveMatrixData( { vector<double>vec1,
   //                        vector<double>vec2,
   //                        vector<double>vec3 } )
   // PLOT.SaveMatrixData( { {1, 3, 5, 6, 1, 4, 5, 4},
   //                        {1, 3, 4, 6, 2, 3, 1, 2},
   //                        {3, 4, 6, 2, 1, 4, 1, 2},
   //                        {2, 5, 5, 4, 6, 7, 4, 4} } )
   std::vector<std::map<std::string, std::string>> MatrixOPT;
   const std::map<std::string, std::string> initial_qt_set{{"term", "qt noraise enhanced font \"Times,22\""},
                                                           {"xlabel", "\'{/Times:Italic=22 x}\'"},
                                                           {"ylabel", "\'{/Times:Italic=22 y}\' rotate by 0"},
                                                           {"tics", "font \"Times,22\""},
                                                           {"grid", "xtics ytics"},
                                                           {"key", "off"},
                                                           {"linestyle 1", " dt 1 lw 2 ps .7 pt 6 lc rgb \"#e41a1c\""},
                                                           {"linestyle 2", " dt 2 lw 2 ps .7 pt 7 lc rgb \"#377eb8\""},
                                                           {"linestyle 3", " dt 4 lw 2 ps .7 pt 8 lc rgb \"#4daf4a\""},
                                                           {"linestyle 4", " dt 5 lw 2 ps .7 pt 9 lc rgb \"#984ea3\""},
                                                           {"linestyle 5", " dt 3 lw 2 ps .7 pt 4 lc rgb \"#ff7f00\""},
                                                           {"linestyle 6", " dt 1 lw 2 ps .7 pt 5 lc rgb '#ED34ED'"},
                                                           {"linestyle 7", " dt 1 lw 2 ps .7 pt 3 lc rgb '#0000FF'"},
                                                           {"style arrow 1", "nohead lc 'red' lw 1 dt 1"},
                                                           {"style arrow 2", "nohead lc 'blue' lw 1 dt 1"},
                                                           {"style arrow 3", "nohead lc \"#4daf4a\" lw 1 dt 1"},
                                                           {"style arrow 4", "nohead lc \"#984ea3\" lw 1 dt 1"},
                                                           {"style arrow 5", "nohead lc \"#ff7f00\" lw 1 dt 1"},
                                                           {"style arrow 6", "nohead lc '#ED34ED' lw 1 dt 1"},
                                                           {"style arrow 7", "nohead lc '#0000FF' lw 1 dt 1"}};

   const std::map<std::string, std::string> def_qt_set{{"term", "qt noraise enhanced font \"Times,22\""},
                                                       {"xlabel", "\'{/Times:Italic=22 x}\'"},
                                                       {"ylabel", "\'{/Times:Italic=22 y}\' rotate by 0"},
                                                       {"tics", "font \"Times,22\""},
                                                       {"grid", "xtics ytics"},
                                                       {"key", "off"},
                                                       {"linestyle 1", " dt 1 lw 2 ps .7 pt 6 lc rgb \"#e41a1c\""},
                                                       {"linestyle 2", " dt 2 lw 2 ps .7 pt 7 lc rgb \"#377eb8\""},
                                                       {"linestyle 3", " dt 4 lw 2 ps .7 pt 8 lc rgb \"#4daf4a\""},
                                                       {"linestyle 4", " dt 5 lw 2 ps .7 pt 9 lc rgb \"#984ea3\""},
                                                       {"linestyle 5", " dt 3 lw 2 ps .7 pt 4 lc rgb \"#ff7f00\""},
                                                       {"linestyle 6", " dt 1 lw 2 ps .7 pt 5 lc rgb '#ED34ED'"},
                                                       {"linestyle 7", " dt 1 lw 2 ps .7 pt 3 lc rgb '#0000FF'"},
                                                       {"style arrow 1", "nohead lc 'red' lw 1 dt 1"},
                                                       {"style arrow 2", "nohead lc 'blue' lw 1 dt 1"},
                                                       {"style arrow 3", "nohead lc \"#4daf4a\" lw 1 dt 1"},
                                                       {"style arrow 4", "nohead lc \"#984ea3\" lw 1 dt 1"},
                                                       {"style arrow 5", "nohead lc \"#ff7f00\" lw 1 dt 1"},
                                                       {"style arrow 6", "nohead lc '#ED34ED' lw 1 dt 1"},
                                                       {"style arrow 7", "nohead lc '#0000FF' lw 1 dt 1"}};

   const std::map<std::string, std::string> initial_eps_set{{"term", "postscript eps enhanced color font \"Helvetica,20\" linewidth 2"},
                                                            {"xlabel", "\'{/Helvetica-Oblique=20 x}\'"},
                                                            {"ylabel", "\'{/Helvetica-Oblique=20 y}\' rotate by 0"},
                                                            {"tics", "font \"Helvetica,20\""},
                                                            {"grid", "xtics ytics"},
                                                            {"key", "off"},
                                                            {"linestyle 1", " dt 1 lw 2 ps .7 pt 6 lc rgb \"#e41a1c\""},
                                                            {"linestyle 2", " dt 1 lw 2 ps .7 pt 7 lc rgb \"#377eb8\""},
                                                            {"linestyle 3", " dt 1 lw 2 ps .7 pt 8 lc rgb \"#4daf4a\""},
                                                            {"linestyle 4", " dt 5 lw 2 ps .7 pt 9 lc rgb \"#984ea3\""},
                                                            {"linestyle 5", " dt 3 lw 2 ps .7 pt 4 lc rgb \"#ff7f00\""},
                                                            {"linestyle 6", " dt 1 lw 2 ps .7 pt 5 lc rgb '#ED34ED'"},
                                                            {"linestyle 7", " dt 1 lw 2 ps .7 pt 3 lc rgb '#0000FF'"},
                                                            {"style arrow 1", "nohead lc 'red' lw 1 dt 1"},
                                                            {"style arrow 2", "nohead lc 'blue' lw 1 dt 1"},
                                                            {"style arrow 3", "nohead lc \"#4daf4a\" lw 1 dt 1"},
                                                            {"style arrow 4", "nohead lc \"#984ea3\" lw 1 dt 1"},
                                                            {"style arrow 5", "nohead lc \"#ff7f00\" lw 1 dt 1"},
                                                            {"style arrow 6", "nohead lc '#ED34ED' lw 1 dt 1"},
                                                            {"style arrow 7", "nohead lc '#0000FF' lw 1 dt 1"}};

   const std::map<std::string, std::string> def_eps_set{{"term", "postscript eps enhanced color font \"Helvetica,20\" linewidth 2"},
                                                        {"xlabel", "\'{/Helvetica-Oblique=20 x}\'"},
                                                        {"ylabel", "\'{/Helvetica-Oblique=20 y}\' rotate by 0"},
                                                        {"tics", "font \"Helvetica,20\""},
                                                        {"grid", "xtics ytics"},
                                                        {"size", "1., 1."},
                                                        {"key", "off"},
                                                        {"linestyle 1", " dt 1 lw 2 ps .7 pt 6 lc rgb \"#e41a1c\""},
                                                        {"linestyle 2", " dt 1 lw 2 ps .7 pt 7 lc rgb \"#377eb8\""},
                                                        {"linestyle 3", " dt 1 lw 2 ps .7 pt 8 lc rgb \"#4daf4a\""},
                                                        {"linestyle 4", " dt 5 lw 2 ps .7 pt 9 lc rgb \"#984ea3\""},
                                                        {"linestyle 5", " dt 3 lw 2 ps .7 pt 4 lc rgb \"#ff7f00\""},
                                                        {"linestyle 6", " dt 1 lw 2 ps .7 pt 5 lc rgb '#ED34ED'"},
                                                        {"linestyle 7", " dt 1 lw 2 ps .7 pt 3 lc rgb '#0000FF'"},
                                                        {"style arrow 1", "nohead lc 'red' lw 1 dt 1"},
                                                        {"style arrow 2", "nohead lc 'blue' lw 1 dt 1"},
                                                        {"style arrow 3", "nohead lc \"#4daf4a\" lw 1 dt 1"},
                                                        {"style arrow 4", "nohead lc \"#984ea3\" lw 1 dt 1"},
                                                        {"style arrow 5", "nohead lc \"#ff7f00\" lw 1 dt 1"},
                                                        {"style arrow 6", "nohead lc '#ED34ED' lw 1 dt 1"},
                                                        {"style arrow 7", "nohead lc '#0000FF' lw 1 dt 1"}};

   const std::vector<std::string> OptList{"using", "u",
                                          "with", "w",
                                          "linecoler", "lc",
                                          "pointtype", "pt",
                                          "pointsize", "ps",
                                          "linewidth", "lw",
                                          "linetype", "lt",
                                          "linestyle", "ls",
                                          "title", "tit",
                                          "arrowstyle", "notitle",
                                          "dashtype", "dt"};

   const std::vector<std::string> OptListWithQuotation{"title", "tit", "notitle"};

   const std::vector<std::string> SetList{"term", "terminal", "angles", "arrow", "autoscale", "bars",
                                          "bmargin", "border", "boxwidth", "cbdata",
                                          "cbdtics", "cblabel", "cbmtics", "cbrange",
                                          "cbtics", "clabel", "clip", "cntrlabel",
                                          "cntrparam", "color", "colorbox", "colorsequence",
                                          "contour", "dashtype", "data", "datafile",
                                          "date_specifiers", "decimalsign", "dgrid3d", "dummy",
                                          "encoding", "errorbars", "fit", "fontpath",
                                          "format", "function", "grid", "hidden3d",
                                          "history", "historysize", "isosamples", "jitter",
                                          "key", "label", "linetype", "link",
                                          "lmargin", "loadpath", "locale", "log",
                                          "logscale", "macros", "mapping", "margin",
                                          "margins", "micro", "minussign", "missing",
                                          "monochrome", "mouse", "mttics", "multiplot",
                                          "mx2tics", "mxtics", "my2tics", "mytics",
                                          "mztics", "nonlinear", "object", "offsets",
                                          "origin", "output", "palette", "parametric", "pm3d",
                                          "polar", "print", "psdir", "raxis",
                                          "rgbmax", "rlabel", "rmargin", "rrange",
                                          "rtics", "samples", "size", "style",
                                          "surface", "table",
                                          "termoption", "tsPlheta", "tics", "ticscale",
                                          "ticslevel", "time_specifiers", "timefmt", "timestamp",
                                          "title", "tmargin", "trange", "ttics",
                                          "urange", "view", "vrange", "x2data",
                                          "x2dtics", "x2label", "x2mtics", "x2range",
                                          "x2tics", "x2zeroaxis", "xdata", "xdtics",
                                          "xlabel", "xmtics", "xrange", "xtics",
                                          "xyplane", "xzeroaxis", "y2data", "y2dtics",
                                          "y2label", "y2mtics", "y2range", "y2tics",
                                          "y2zeroaxis", "ydata", "ydtics", "ylabel",
                                          "ymtics", "yrange", "ytics", "yzeroaxis",
                                          "zdata", "zdtics", "zero", "zeroaxis",
                                          "zlabel", "zmtics", "zrange", "ztics", "zzeroaxis",
                                          "linestyle 1", "linestyle 2", "linestyle 3", "linestyle 4", "linestyle 5", "linestyle 6",
                                          "style arrow 1", "style arrow 2", "style arrow 3", "style arrow 4", "style arrow 5", "style arrow 6",
                                          "style arrow 7", "style arrow 8"};

   template <class T>
   std::string rgb(T r, T g, T b) {
      std::stringstream ss;
      ss << "rgb \'#" << std::hex << ((int)r << 16 | (int)g << 8 | (int)b) << "\'";
      return ss.str();
   }

   std::string rgb(const std::vector<double>& abc) {
      return rgb(abc[0], abc[1], abc[2]);
   }

   void manageSet(std::map<std::string, std::string> mp, bool set = true) {
      std::string prefix = set ? "set " : "unset ";
      for (const auto& list : SetList) {
         auto it = mp.find(list);
         if (it != mp.end()) {
            std::string str(prefix + it->first + " " + it->second);
            std::fprintf(proc, "%s\n", str.c_str());
   #ifdef debug_GNUPLOT
            std::cout << str << std::endl;
   #endif
         }
      }
   };

   void Set(std::map<std::string, std::string> mp) {
      manageSet(mp, true);
   };

   void UnSet(std::map<std::string, std::string> mp) {
      manageSet(mp, false);
   };

   void InitGnuplot(const std::map<std::string, std::string>& mp = {}) {
   #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__) || defined(_WIN64)
      proc = _popen("gnuplot -persist", "w");
   #else
      proc = popen("gnuplot -persist", "w");
   #endif

      if (!proc)
         std::cout << "ERROR" << __func__ << std::endl;
      Set(initial_qt_set);
      if (!mp.empty()) Set(mp);

      std::string str = "\nMATLAB = \"defined (0  0.0 0.0 0.5, 1  0.0 0.0 1.0, 2  0.0 0.5 1.0, 3  0.0 1.0 1.0, 4  0.5 g1.0 0.5, 5  1.0 1.0 0.0, 6  1.0 0.5 0.0, 7  1.0 0.0 0.0, 8  0.5 0.0 0.0 )\"";
      Set({{"macros", str}});
   };

   GNUPLOT() : counter(0), nplot(0), n_decld_plot(0) {
      InitGnuplot();
   };

   GNUPLOT(const std::map<std::string, std::string>& mp)
       : counter(0), nplot(0), n_decld_plot(0) {
      InitGnuplot(mp);
   };

   ~GNUPLOT() {
      std::fprintf(proc, "set term qt 0 close\n");
   #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__) || defined(_WIN64)
      _pclose(proc);
   #else
      pclose(proc);
   #endif
      std::cout << "GNUPLOT is destructed" << std::endl;
   };

   void Clear() {
      Data.clear();
      OPT.clear();
      VectorData.clear();
      VectorOPT.clear();
      MatrixData.clear();
      MatrixOPT.clear();
      SplotData.clear();
      SplotOPT.clear();
   };
   //============================================
   std::string map_to_string(const std::map<std::string, std::string>& mp) {
      std::string opt;
      bool found = false;
      for (const auto& list : OptList) {
         if (mp.find(list) != mp.end()) {

            for (const auto& needquotation : OptListWithQuotation) {
               if (mp.find(list)->first == needquotation) {
                  found = true;
                  break;
               }
            }

            if (found) {
               opt = opt + " " + mp.find(list)->first + " \'" + mp.find(list)->second + "\'";
            } else {
               opt = opt + " " + mp.find(list)->first + " " + mp.find(list)->second;
            }
         }
      }
      return opt;
   };

   void PlotAfterDecl(std::vector<std::vector<std::vector<double>>> v) {
      for (const auto& u : v) {
         for (size_t i = 0; i < u[0].size(); i++)
            std::fprintf(proc, "%f  %f\n", u[0][i], u[1][i]);
         std::fprintf(proc, "%s\n", "e");
      }
      std::fflush(proc);
   };

   void GetDeclString(int s, std::vector<std::string>& opt, std::string& str) {
      for (auto i = 0; i < s; i++) {
         str = str + ((std::string)(nplot > 0 ? "," : "plot")) + " " + "'-'" + " " + opt[(int)opt.size() - i % (int)opt.size() - 1];
         nplot++;
      }
   };
   ////////////////////////////////////////////////////
   // ex.
   // PLOT.SaveData({vector<double>x})
   // PLOT.SaveData({vector<double>x, vector<double>y})
   // PLOT.SaveData({vector<double>x, vector<double>y, vector<double>z})
   void SaveData(const std::vector<double>& v) {
      Data.push_back({v});
      OPT.push_back(std::map<std::string, std::string>{{"ls", std::to_string((int)(Data.size() % 7 + 1))}});
   };
   void SaveData(const std::vector<double>& v, const std::map<std::string, std::string>& mp) {
      Data.push_back({v});
      OPT.push_back(mp);
   };
   void SaveData(const std::vector<std::vector<double>>& v) {
      Data.push_back(v);
      OPT.push_back(std::map<std::string, std::string>{{"ls", std::to_string((int)(Data.size() % 7 + 1))}});
   };
   void SaveData(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
      auto w = v;
      if (mp.find("loop") != mp.end()) {
         w.push_back(v[0]);
         Data.push_back(w);
      } else {
         Data.push_back(v);
      }
      OPT.push_back(mp);
   };
   ////////////////////////////////////////////////////
   // ex.
   // PLOT.SaveSplotData( { {1, 1, 2}, {2, 1, 2}, {3, 1, 2}, {4, 1, 1} },
   //                       {1, 2, 2}, {2, 2, 2}, {3, 2, 2}, {4, 2, 2} },
   //                       {1, 3, 2}, {2, 3, 2}, {3, 3, 2}, {4, 3, 1} } })
   void SaveSplotData(const std::vector<std::vector<std::vector<double>>>& v) {
      SplotData.push_back(v);
      SplotOPT.push_back(std::map<std::string, std::string>{{"ls", std::to_string((int)(Data.size() % 5 + 1))}});
   };
   void SaveSplotData(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
      SplotData.push_back(v);
      SplotOPT.push_back(mp);
   };
   ////////////////////////////////////////////////////
   // ex.
   // PLOT.SaveVectorData({{vector<double>x1, vector<double>y1}, {vector<double>x2, vector<double>y2}})
   // PLOT.SaveVectorData( { {3, 4, 2}, {3, 4, 2} },
   //                        {3, 4, 2}, {3, 4, 2} },
   //                        {3, 4, 2}, {3, 4, 2} } })
   void SaveVectorData(const std::vector<std::vector<std::vector<double>>>& v) {
      VectorData.push_back(v);
      VectorOPT.push_back(std::map<std::string, std::string>{{"ls", std::to_string((int)(Data.size() % 5 + 1))}});
   };
   void SaveVectorData(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
      VectorData.push_back(v);
      VectorOPT.push_back(mp);
   };
   ////////////////////////////////////////////////////
   // ex.
   // PLOT.SaveMatrixData( { vector<double>vec1,
   //                        vector<double>vec2,
   //                        vector<double>vec3 } )
   // PLOT.SaveMatrixData( { {1, 3, 5, 6, 1, 4, 5, 4},
   //                        {1, 3, 4, 6, 2, 3, 1, 2},
   //                        {3, 4, 6, 2, 1, 4, 1, 2},
   //                        {2, 5, 5, 4, 6, 7, 4, 4} } )
   void SaveMatrixData(const std::vector<std::vector<double>>& v) {
      MatrixData.push_back(v);
      MatrixOPT.push_back(std::map<std::string, std::string>{{"ls", std::to_string((int)(Data.size() % 5 + 1))}});
   };
   void SaveMatrixData(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
      MatrixData.push_back(v);
      MatrixOPT.push_back(mp);
   };
   ///////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////
   void declare_Plot1D(const int n) {
      if (n_decld_plot == 0)
         decld_str = "plot '-' " + map_to_string(OPT[n]);
      else
         decld_str = decld_str + ", '-' " + map_to_string(OPT[n]);

      n_decld_plot++;
   };
   void declare_Plot2D(const int n) {
      declare_Plot1D(n);
   };
   void declare_Plot3D(const int n) {
      if (n_decld_plot == 0)
         decld_str = "splot '-' " + map_to_string(OPT[n]);
      else
         decld_str = decld_str + ", '-' " + map_to_string(OPT[n]);

      n_decld_plot++;
   };
   //////////////
   void declare_Splot(const int n) {
      if (n_decld_plot == 0)
         decld_str = "splot '-' " + map_to_string(SplotOPT[n]);
      else
         decld_str = decld_str + ", '-' " + map_to_string(SplotOPT[n]);

      n_decld_plot++;
   };
   //////////////
   void declare_VectorPlot(const int n) {
      if (n_decld_plot == 0)
         decld_str = "plot '-' w vec " + map_to_string(VectorOPT[n]);
      else
         decld_str = decld_str + ", '-' w vec " + map_to_string(VectorOPT[n]);

      n_decld_plot++;
   };
   void declare_VectorPlot3D(const int n) {
      if (n_decld_plot == 0)
         decld_str = "splot '-' w vec " + map_to_string(VectorOPT[n]);
      else
         decld_str = decld_str + ", '-' w vec " + map_to_string(VectorOPT[n]);

      n_decld_plot++;
   };
   //////////////
   void declare_MatrixPlot(const int n) {
      if (n_decld_plot == 0)
         decld_str = "plot '-' matrix " + map_to_string(MatrixOPT[n]);
      else
         decld_str = decld_str + ", '-' matrix" + map_to_string(MatrixOPT[n]);

      n_decld_plot++;
   };
   //////////////
   void declare_Plot1D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_Plot1D(n);
   };
   void declare_Plot2D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_Plot2D(n);
   };
   void declare_Plot3D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_Plot3D(n);
   };
   void declare_Splot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_Splot(n);
   };
   void declare_VectorPlot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_VectorPlot(n);
   };
   void declare_VectorPlot3D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_VectorPlot3D(n);
   };
   void declare_MatrixPlot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         declare_MatrixPlot(n);
   };
   //////////////
   void declare_Plot1D() {
      for (size_t n = 0; n < Data.size(); n++)
         declare_Plot1D(n);
   };
   void declare_Plot2D() {
      for (size_t n = 0; n < Data.size(); n++)
         declare_Plot2D(n);
   };
   void declare_Plot3D() {
      for (size_t n = 0; n < Data.size(); n++)
         declare_Plot3D(n);
   };
   void declare_Splot() {
      for (size_t n = 0; n < Data.size(); n++)
         declare_Splot(n);
   };
   void declare_VectorPlot() {
      for (size_t n = 0; n < VectorData.size(); n++)
         declare_VectorPlot(n);
   };
   void declare_VectorPlot3D() {
      for (size_t n = 0; n < VectorData.size(); n++)
         declare_VectorPlot3D(n);
   };
   void declare_MatrixPlot() {
      for (size_t n = 0; n < Data.size(); n++)
         declare_MatrixPlot(n);
   };
   ///////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////
   void proc_Plot1D(const int n) {
      for (size_t i = 0; i < Data[n].size(); i++)
         std::fprintf(proc, "%f\n", Data[n][i][0]);
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_Plot1D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_Plot1D(n);
   };  // plot multiple data set
   void Plot1D(const int n) {
      declare_Plot1D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot1D(n);
      n_decld_plot = 0;
   };
   void Plot1D(const std::vector<double>& v) {
      SaveData({v});
      int n((int)Data.size() - 1);
      declare_Plot1D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot1D(n);
      n_decld_plot = 0;
   };
   void Plot1D(const std::vector<double>& v, const std::map<std::string, std::string>& mp) {
      SaveData({v}, mp);
      int n((int)Data.size() - 1);
      declare_Plot1D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot1D(n);
      n_decld_plot = 0;
   };
   void Plot1D() {
      std::vector<int> N(Data.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      declare_Plot1D(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot1D(N);
      n_decld_plot = 0;
   };
   //////////////////////////////////////////////////////
   void proc_Plot2D(const int n) {
      for (size_t i = 0; i < Data[n].size(); i++)
         std::fprintf(proc, "%f %f\n", Data[n][i][0], Data[n][i][1]);
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_Plot2D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_Plot2D(n);
   };  // plot multiple data set
   void Plot2D(const int n) {
      declare_Plot2D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot2D(n);
      n_decld_plot = 0;
   };
   void Plot2D(const std::vector<std::vector<double>>& v) {
      SaveData(v);
      int n((int)Data.size() - 1);
      declare_Plot2D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot2D(n);
      n_decld_plot = 0;
   };
   void Plot2D(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
      SaveData(v, mp);
      int n((int)Data.size() - 1);
      declare_Plot2D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot2D(n);
      n_decld_plot = 0;
   };
   void Plot2D() {
      std::vector<int> N(Data.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      declare_Plot2D(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot2D(N);
      n_decld_plot = 0;
   };
   void plot2d() {
      Plot2D();
   };
   //////////////////////////////////////////////////////
   void proc_Plot3D(const int n) {
      for (size_t i = 0; i < Data[n].size(); i++)
         std::fprintf(proc, "%f %f %f\n", Data[n][i][0], Data[n][i][1], Data[n][i][2]);
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_Plot3D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_Plot3D(n);
   };  // plot multiple data set
   void Plot3D(const int n) {
      declare_Plot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(n);
      n_decld_plot = 0;
   };
   void Plot3D(const std::vector<std::vector<double>>& v) {
      SaveData(v);
      int n((int)Data.size() - 1);
      declare_Plot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(n);
      n_decld_plot = 0;
   };
   void Plot3D(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
      SaveData(v, mp);
      int n((int)Data.size() - 1);
      declare_Plot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(n);
      n_decld_plot = 0;
   };
   void Plot3D() {
      std::vector<int> N(Data.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      declare_Plot3D(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(N);
      n_decld_plot = 0;
   };
   //////////////////////////////////////////////////////
   void proc_Splot(const int n) {
      for (size_t i = 0; i < SplotData[n].size(); i++) {
         for (size_t j = 0; j < SplotData[n][i].size(); j++) {
            switch ((int)SplotData[n][i][j].size()) {
               case 3:
                  std::fprintf(proc, "%f %f %f\n", SplotData[n][i][j][0], SplotData[n][i][j][1], SplotData[n][i][j][2]);
                  break;
               case 4:
                  std::fprintf(proc, "%f %f %f %f\n", SplotData[n][i][j][0], SplotData[n][i][j][1], SplotData[n][i][j][2], SplotData[n][i][j][3]);
                  break;
               case 5:
                  std::fprintf(proc, "%f %f %f %f %f\n", SplotData[n][i][j][0], SplotData[n][i][j][1], SplotData[n][i][j][2], SplotData[n][i][j][3], SplotData[n][i][j][4]);
                  break;
               case 6:
                  std::fprintf(proc, "%f %f %f %f %f %f\n", SplotData[n][i][j][0], SplotData[n][i][j][1], SplotData[n][i][j][2], SplotData[n][i][j][3], SplotData[n][i][j][4], SplotData[n][i][j][5]);
                  break;
            }
         }
         std::fprintf(proc, " \n");
      }
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // splot
   void proc_Splot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_Splot(n);
   };  // plot multiple data set
   void Splot(const int n) {
      declare_Splot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Splot(n);
      n_decld_plot = 0;
   };
   void Splot(const std::vector<std::vector<std::vector<double>>>& v) {
      SaveSplotData(v);
      int n((int)SplotData.size() - 1);
      declare_Splot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Splot(n);
      n_decld_plot = 0;
   };
   void Splot(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
      SaveSplotData(v, mp);
      int n((int)SplotData.size() - 1);
      declare_Splot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Splot(n);
      n_decld_plot = 0;
   };
   void Splot() {
      std::vector<int> N(SplotData.size());
      for (size_t n = 0; n < SplotData.size(); n++)
         N[n] = n;
      declare_Splot(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Splot(N);
      n_decld_plot = 0;
   };
   //////////////////////////////////////////////////////
   // plot vector in 2d
   void proc_VectorPlot(const int n) {
      if (VectorOPT[n]["lc"] == "palette") {
         for (size_t i = 0; i < VectorData[n].size(); i++)
            std::fprintf(proc, "%f %f %f %f %f\n",
                         VectorData[n][i][0][0], VectorData[n][i][0][1], /*座標x*/
                         VectorData[n][i][1][0], VectorData[n][i][1][1], /*ベクトルdx*/
                         VectorData[n][i][2][0]                          /*palette value*/
            );
      } else {
         for (size_t i = 0; i < VectorData[n].size(); i++)
            std::fprintf(proc, "%f %f %f %f\n",
                         VectorData[n][i][0][0], VectorData[n][i][0][1], /*座標x*/
                         VectorData[n][i][1][0], VectorData[n][i][1][1]  /*ベクトルdx*/
            );
      }
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_VectorPlot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_VectorPlot(n);
   };  // plot multiple data set
   void VectorPlot(const int n) {
      declare_VectorPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot(n);
      n_decld_plot = 0;
   };
   void VectorPlot(const std::vector<std::vector<std::vector<double>>>& v) {
      SaveVectorData(v);
      int n((int)VectorData.size() - 1);
      declare_VectorPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot(n);
      n_decld_plot = 0;
   };
   void VectorPlot(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
      SaveVectorData(v, mp);
      int n((int)VectorData.size() - 1);
      declare_VectorPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot(n);
      n_decld_plot = 0;
   };
   void VectorPlot() {
      std::vector<int> N(VectorData.size());
      for (size_t n = 0; n < VectorData.size(); n++)
         N[n] = n;
      declare_VectorPlot(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot(N);
      n_decld_plot = 0;
   };
   //////////////////////////////////////////////////////
   // plot vector in 3d
   void proc_VectorPlot3D(const int n) {
      if (VectorOPT[n]["lc"] == "palette") {
         for (size_t i = 0; i < VectorData[n].size(); i++)
            std::fprintf(proc, "%f %f %f %f %f %f %f\n",
                         VectorData[n][i][0][0], VectorData[n][i][0][1], VectorData[n][i][0][2], /*座標x*/
                         VectorData[n][i][1][0], VectorData[n][i][1][1], VectorData[n][i][1][2], /*ベクトルdx*/
                         VectorData[n][i][2][0]                                                  /*palette value*/
            );
      } else {
         for (size_t i = 0; i < VectorData[n].size(); i++)
            std::fprintf(proc, "%f %f %f %f %f %f\n",
                         VectorData[n][i][0][0], VectorData[n][i][0][1], VectorData[n][i][0][2], /*座標x*/
                         VectorData[n][i][1][0], VectorData[n][i][1][1], VectorData[n][i][1][2]  /*ベクトルdx*/
            );
      }
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_VectorPlot3D(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_VectorPlot3D(n);
   };  // plot multiple data set
   void VectorPlot3D(const int n) {
      declare_VectorPlot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot3D(n);
      n_decld_plot = 0;
   };
   void VectorPlot3D(const std::vector<std::vector<std::vector<double>>>& v) {
      SaveVectorData(v);
      int n((int)VectorData.size() - 1);
      declare_VectorPlot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot3D(n);
      n_decld_plot = 0;
   };
   void VectorPlot3D(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
      SaveVectorData(v, mp);
      int n((int)VectorData.size() - 1);
      declare_VectorPlot3D(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot3D(n);
      n_decld_plot = 0;
   };
   void VectorPlot3D() {
      std::vector<int> N(VectorData.size());
      for (size_t n = 0; n < VectorData.size(); n++)
         N[n] = n;
      declare_VectorPlot3D(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_VectorPlot3D(N);
      n_decld_plot = 0;
   };
   //////////////////////////////////////////////////////
   void Plot_All() {
      std::vector<int> N(Data.size());
      std::vector<int> V(VectorData.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      for (size_t n = 0; n < VectorData.size(); n++)
         V[n] = n;
      declare_Plot2D(N);
      declare_VectorPlot(V);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot2D(N);
      proc_VectorPlot(V);
      n_decld_plot = 0;
   };
   void Plot3D_All() {
      n_decld_plot = 0;  // if 0 is set, plot will be initialized
      std::vector<int> N(Data.size());
      std::vector<int> S(SplotData.size());
      std::vector<int> V(VectorData.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      for (size_t n = 0; n < SplotData.size(); n++)
         S[n] = n;
      for (size_t n = 0; n < VectorData.size(); n++)
         V[n] = n;
      declare_Plot3D(N);
      declare_Splot(S);
      declare_VectorPlot3D(V);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(N);
      proc_Splot(S);
      proc_VectorPlot3D(V);
   };
   void plot3d() {       // 初めからプロットし直す
      n_decld_plot = 0;  // if 0 is set, plot will be initialized: start proc as 'plot'
      std::vector<int> N(Data.size());
      std::vector<int> S(SplotData.size());
      std::vector<int> V(VectorData.size());
      for (size_t n = 0; n < Data.size(); n++)
         N[n] = n;
      for (size_t n = 0; n < SplotData.size(); n++)
         S[n] = n;
      for (size_t n = 0; n < VectorData.size(); n++)
         V[n] = n;
      declare_Plot3D(N);  // declaring data
      declare_Splot(S);
      declare_VectorPlot3D(V);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_Plot3D(N);  // outputting stored data
      proc_Splot(S);
      proc_VectorPlot3D(V);
   };
   //////////////////////////////////////////////////////
   void proc_MatrixPlot(const int n) {
      for (size_t i = 0; i < MatrixData[n].size(); i++) {
         for (size_t j = 0; j < MatrixData[n][i].size(); j++)
            std::fprintf(proc, "%f ", MatrixData[n][i][j]);
         std::fprintf(proc, "\n ");
      }
      std::fprintf(proc, "%s\n", "e");
      std::fflush(proc);
   };  // plot single data
   void proc_MatrixPlot(const std::vector<int>& N) {
      for (size_t n = 0; n < N.size(); n++)
         proc_MatrixPlot(n);
   };  // plot multiple data set
   void MatrixPlot(const int n) {
      declare_MatrixPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_MatrixPlot(n);
      n_decld_plot = 0;
   };
   void MatrixPlot(const std::vector<std::vector<double>>& v) {
      SaveMatrixData(v);
      int n((int)MatrixData.size() - 1);
      declare_MatrixPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_MatrixPlot(n);
      n_decld_plot = 0;
   };
   void MatrixPlot(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
      SaveMatrixData(v, mp);
      int n((int)MatrixData.size() - 1);
      declare_MatrixPlot(n);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_MatrixPlot(n);
      n_decld_plot = 0;
   };
   void MatrixPlot() {
      std::vector<int> N(MatrixData.size());
      for (size_t n = 0; n < MatrixData.size(); n++)
         N[n] = n;
      declare_MatrixPlot(N);
      std::fprintf(proc, "%s\n", decld_str.c_str());
      proc_MatrixPlot(N);
      n_decld_plot = 0;
   };
};

//////////////////////////////////////////////////////

// void Plot3DAll(){
//   std::vector<int> N(Data.size());
//   for(auto n=0; n<Data.size(); n++)
//     N[n] = n;
//   std::vector<int> M(VectorData.size());
//   for(auto n=0; n<VectorData.size(); n++)
//     M[n] = n;

//   DeclarePlot3D(N);
//   DeclareVectorPlot3D(M);

// };

//////////////////////////////////////////////////////
// ListPlot function

GNUPLOT* Plot2D(const std::vector<std::vector<double>>& v) {
   GNUPLOT* plot = new GNUPLOT;
   plot->Plot2D(v);
   std::cin.ignore();
   return plot;
};

GNUPLOT* Plot3D(const std::vector<std::vector<double>>& v) {
   GNUPLOT* plot = new GNUPLOT;
   plot->Plot3D(v);
   std::cin.ignore();
   return plot;
};

GNUPLOT* Plot3D(const std::vector<std::vector<double>>& v, const std::map<std::string, std::string>& mp) {
   GNUPLOT* plot = new GNUPLOT;
   plot->Plot3D(v, mp);
   std::cin.ignore();
   return plot;
};

GNUPLOT* VectorPlot3D(const std::vector<std::vector<std::vector<double>>>& v) {
   GNUPLOT* plot = new GNUPLOT;
   plot->VectorPlot3D(v);
   std::cin.ignore();
   return plot;
};

GNUPLOT* VectorPlot3D(const std::vector<std::vector<std::vector<double>>>& v, const std::map<std::string, std::string>& mp) {
   GNUPLOT* plot = new GNUPLOT;
   plot->VectorPlot3D(v, mp);
   std::cin.ignore();
   return plot;
};

#endif

////////////////////
///// EXAMPLES /////
////////////////////

/** std::vector<std::vector<double>> data;
 * double x, y;
 * int N=100;
 * for(int i=0; i<N; i++)
 *   {
 *     x = (double)i -50;
 *     y = pow(cosh(x/(2.*pi)),-2);
 *     std::vector<double> X{x,y};
 *     data.push_back(X);
 *   }
 * cout << data << endl;
 * GNUPLOT PLOT;
 * PLOT.Plot2D(data,{{"ls","1"},{"w","lp"}});
 */

/** std::vector<std::vector<std::vector<double>>> data;
 * double x, y, z;
 * int N=20;
 * for(int i=0; i<N; i++)
 *   {
 *     std::vector<std::vector<double>> subdata;
 *     x = (double)i - N/2;
 *     for(int j=0; j<N; j++)
 * 	{
 * 	  y = (double)j - N/2;
 * 	  z = pow(cosh((x+y)/(2.*pi)),-2);
 * 	  std::vector<double> X{x,y,z};
 * 	  subdata.push_back(X);
 * 	}
 *     data.push_back(subdata);
 *   }
 * cout << data << endl;
 * GNUPLOT PLOT;
 * PLOT.Set({{"hidden3d",""}});
 * PLOT.Splot(data,{{"w","l"}});
 */
////////////

/** std::vector<std::vector<std::vector<double>>> data;
 * double x, y;
 * int N=20;
 * for(int i=0; i<N; i++)
 *   for(int j=0; j<N; j++)
 *     {
 * 	x = 3.*i/N - 1.5;
 * 	y = 3.*j/N - 1.5;
 * 	std::vector<double> X{x,y};
 * 	std::vector<double> dX{y, x - pow(x,3)};
 * 	data.push_back({ X, dX/5. , {Norm(dX)} });
 *     }
 * GNUPLOT PLOT;
 * PLOT.VectorPlot(data ,{{"lc","palette"}});
 */
///////////////////

/** std::vector<std::vector<std::vector<double>>> data;
 * double x, y, z;
 * int N=10;
 * for(int i=0; i<N; i++)
 *   for(int j=0; j<N; j++)
 *     for(int k=0; k<N; k++)
 * 	{
 * 	  x = 3.*i/N - 1.5;
 * 	  y = 3.*j/N - 1.5;
 * 	  z = 3.*k/N - 1.5;
 * 	  std::vector<double> X{x,y,z};
 * 	  std::vector<double> dX{y, x - pow(x,3), z - pow(z,3)};
 * 	  data.push_back({ X, dX/5. , {Norm(dX)} });
 * 	}
 * GNUPLOT PLOT;
 * PLOT.VectorPlot3D(data ,{{"lc","palette"}});
 */

///////////////////

/** vector<vector<double>> data(100,vector<double>(100,0.));
 * for(auto i=0; i<100; i++)
 *   for(auto j=0; j<100; j++)
 *     {
 *       double x = pi/99 * i;
 *       double y = pi/99 * j;
 *       data[i][j] = std::sin(pow(y,2) + x);
 *     }
 * GNUPLOT mat;
 * //mat.Set({{"palette","model RGB rgbformulae 33,13,10"},{"xrange","[0:100]"},{"yrange","[0:100]"}});
 * mat.Set({{"xrange","[0:100]"},{"yrange","[0:100]"}});
 * mat.SaveMatrixData(data,{{"w","image"}});
 * mat.MatrixPlot();
 */

/**  ex.)
 * 20200111
 * IMPORT BATHYMETRY DATA. THEN SPLOT WITH LABELD POINTS.
 */

/** std::vector<std::vector<double>> topo=Import("./data.asc");
 * std::vector<std::vector<double>> z=Take(topo,{0,(int)topo.size()-1,15},{0,(int)topo[0].size()-1,15});
 * std::vector<std::vector<double>> x(z), y(z);
 * std::vector<std::vector<std::vector<double>>> data3d, data3d2;
 * for(auto i = 0; i<x.size(); i++)
 *   {
 *     std::vector<std::vector<double>> data2d, data2d2;
 *     for(auto j = 0; j<x[i].size(); j++)
 *       {
 * 	x[i][j] = (double)100.*(2.*j/(x[i].size()-1.) - 1.);
 * 	y[i][j] = -(double)100.*(2.*i/(x.size()-1.) - 1.);
 * 	data2d.push_back({x[i][j],y[i][j]});
 * 	data2d2.push_back({x[i][j],y[i][j],z[i][j]});
 *       }
 *     data3d.push_back(data2d);
 *     data3d2.push_back(data2d2);
 *   }
 * RBF_multiquadric mult(0.);
 * RBF_interp interp(Flatten(data3d),Flatten(z),mult,0);
 * GNUPLOT PLOT;
 * PLOT.Set({{"contour","base"}});
 * std::vector<std::vector<std::vector<double>>> bot;
 * int n = 0;
 * std::cout << "n = " << n << std::endl;
 * for(auto j=0; j<U.B[n].node_size_row; j++){
 *   std::vector<std::vector<double>> datavec;
 *   for(auto i=0; i<U.B[n].node_size_col; i++){
 *     double x1 = (double)100.*(2.*i/(U.B[n].node_size_col - 1.) - 1.);
 *     double x2 = (double)100.*(2.*j/(U.B[n].node_size_row - 1.) - 1.);
 *     datavec.push_back({x1,x2,interp({x1,x2})});
 *   }
 *   bot.push_back(datavec);
 *  }
 * U.set_init_value(n,bot);
 * PLOT.SaveSplotData(bot,{{"w","l"}});
 * PLOT.SaveSplotData(bot,{{"using","1:2:3:(sprintf(\"(\%d, \%d)\", \$1, \$2))"},{"w","labels point"}});
 * PLOT.Splot();
 */
///////////////////////////

////////////////////////////////////////////
/* EXAMPLE WITH BOUNDARY ELEMENT METHOD */
///////////////////////////////////////////
/** for(auto m=0; m < U.B.size(); m++){
 *   for(auto n=0; n < U.B[0].el.size(); n++){
 * 	std::vector<std::vector<double>> X(U.B[m].el[n].xyz.size(),std::vector<double>(U.B[m].el[n].xyz[0].size()));
 * 	std::vector<std::vector<double>> Y(U.B[m].el[n].xyz.size(),std::vector<double>(U.B[m].el[n].xyz[0].size()));
 * 	std::vector<std::vector<double>> Z(U.B[m].el[n].xyz.size(),std::vector<double>(U.B[m].el[n].xyz[0].size()));
 * 	for(auto i = 0; i<U.B[m].el[n].xyz.size(); i++){
 * 	  for(auto j = 0; j<U.B[m].el[n].xyz[i].size(); j++){
 * 	    X[i][j] = U.B[m].el[n].xyz[i][j][0];
 * 	    Y[i][j] = U.B[m].el[n].xyz[i][j][1];
 * 	    Z[i][j] = U.B[m].el[n].xyz[i][j][2];
 * 	  }
 * 	}
 * 	ParametricInterpolation interpX(X, 3);
 * 	ParametricInterpolation interpY(Y, 3);
 * 	ParametricInterpolation interpZ(Z, 3);
 * 	std::vector<std::vector<std::vector<double>>> mat;
 * 	int N(10);
 * 	l=0;
 * 	for(auto i = -N; i<=N; i++){
 * 	  std::vector<std::vector<double>> vec;
 * 	  for(auto j = -N; j<=N; j++){
 * 	    vec.push_back({interpX((double)i/N,(double)j/N),
 * 			   interpY((double)i/N,(double)j/N),
 * 			   interpZ((double)i/N,(double)j/N),
 * 			   (double)l++});
 * 	  }
 * 	  mat.push_back(vec);
 * 	}
 * 	//	PLOT5.SaveSplotData(mat,{{"using","1:2:3:(sprintf(\"\%i\", $4))"},{"w","labels point"}});
 * 	PLOT5.SaveSplotData(mat,{{"w","l"}});
 * 	PLOT5.Splot();
 * 	std::cin.ignore();
 *   }
 * }
 * PLOT5.Clear();
 */
