#pragma once

#include <string>
#include <unordered_map>
#include <iostream>
#include <sstream>



class Parser
{
    struct val_desc
    {
      std::string value;
      std::string description;
      bool exists;
    };
    bool parsed_ = false;

    std::unordered_map<std::string, val_desc> values_;

    int process_help_notfound (std::string opt)
      {
        bool is_help = (opt.find ("help") != opt.npos);
        if (is_help || !values_.count (opt))
          {
            if (!is_help)
                std::cerr << "ERROR: illegal option detected: " << opt << std::endl;
            std::cout << "Available options:" << std::endl;
            for (auto &&X : values_)
              {
                auto &Name = X.first;
                auto &Desc = X.second;
                std::cout << ("-" + Name) << " ( = " << Desc.value
                          << ") : " << Desc.description << std::endl;
              }
            return -1;
          }
        values_[opt].exists = true;
        return 0;
      }

public:
  Parser() = default;

  void print_values (void)
    {
      for (auto &&X : values_)
        {
          auto &name = X.first;
          auto &desc = X.second;
          std::cout << ("-" + name) << " = " << desc.value  << std::endl;  
        }
    }

  template <typename T>
  void add(std::string name, T defval, std::string description = "")
    {
      std::ostringstream os;
      os << defval;
      auto &desc = values_[name];
      desc.value = os.str ();
      desc.description = description;
      desc.exists = true;
    }

  template <typename T> int get(std::string name, T& res) const
    {
      if (!parsed_)
        {
          std::cerr << "ERROR: commaned line is not parsed yet" << std::endl;
          return -1;
        }
      if (!values_.count (name))
        {
          std::cerr << "ERROR: option " << name << " is not present" << std::endl;
          return -1;
        }
      const auto &desc = values_.find (name)->second;
      const auto &val = desc.value;
      std::istringstream is {val};
      is >> res;
      return 0;
    }

  bool exists(std::string name) const
    {
      if (!parsed_)
        return false;
      if (!values_.count (name))
        return false;
      const auto &desc = values_.find (name)->second;
      return desc.exists;
    }

  int parse(int argc, char **argv)
    {
      int res = 0;
      for (int i = 1; i < argc; ++i) 
        {
          std::string_view opt_view = argv[i];
          auto trim_pos = opt_view.find_first_not_of ('-');
          opt_view.remove_prefix (trim_pos);
          trim_pos = opt_view.find ('=');
          bool have_val = (trim_pos != opt_view.npos);
          std::string_view val_view = opt_view;
          if (have_val)
            opt_view.remove_suffix (opt_view.size () - trim_pos);
          std::string opt{opt_view};
          res = process_help_notfound (opt);
          if  (res < 0)
            return res;
          if (have_val)
            {
              val_view.remove_prefix (trim_pos + 1);
              values_[opt].value = val_view;
            }
        }
      parsed_ = true;
      return 0;
    }

  bool parsed() const noexcept { return parsed_; }
};
