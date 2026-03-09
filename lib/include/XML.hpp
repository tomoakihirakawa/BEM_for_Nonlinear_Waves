#ifndef XML_H
#define XML_H

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

struct XMLElement {
   std::string name;
   std::map<std::string, std::string> attributes;
   std::string content;
   std::function<std::string()> content_func;
   std::function<void(std::stringstream &ss)> writer;
   std::vector<XMLElement *> elements;
   std::vector<std::shared_ptr<XMLElement>> elements_shared;
   /* ------------------------------------------------------ */
   XMLElement(const std::string_view &nameIN) : name(nameIN){};
   XMLElement(const std::string_view &nameIN, const std::map<std::string, std::string> &attributesIN) : name(nameIN), attributes(attributesIN){};
   ~XMLElement(){
       // std::cout << this->name << " is deleted" << std::endl;
   };
   /* ------------------------------------------------------ */
   // std::string getStartTag() const {
   //    std::stringstream ss;
   //    for (const auto &[key, value] : attributes)
   //       ss << " " << key << "=\"" << value << "\"";
   //    return "<" + name + ss.str() + ">";
   // }

   std::string getStartTag() const {
      std::string tag = "<" + name;
      for (const auto &[key, value] : attributes)
         tag += " " + key + "=\"" + value + "\"";
      tag += ">";
      return tag;
   }

   std::string getEndTag() const { return "</" + name + ">"; }
   std::tuple<std::string, std::string> getTags() const { return {this->getStartTag(), this->getEndTag()}; }
   void add(XMLElement *elem) { this->elements.emplace_back(elem); }
   void add(std::shared_ptr<XMLElement> elem) { this->elements_shared.emplace_back(elem); }

   template <typename T>
   T &write(T &ofs) const {
      if (this->elements.empty() && this->elements_shared.empty()) {
         if (writer) {
            std::stringstream ss;
            writer(ss);
            ofs << this->getStartTag() << "\n" + ss.str() << "\n"
                << this->getEndTag() << "\n";
         } else
            ofs << this->getStartTag() << "\n"
                << this->getEndTag() << "\n";
      } else {
         ofs << this->getStartTag() << "\n";
         for (const auto &elem : this->elements)
            elem->write(ofs);
         for (const auto &elem : this->elements_shared)
            elem->write(ofs);
         ofs << this->getEndTag() << "\n";
      }
      return ofs;
   }
};

#endif