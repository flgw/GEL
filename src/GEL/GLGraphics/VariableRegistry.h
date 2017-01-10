//
// Created by flgw on 1/10/17.
//

#ifndef GEL_VARIABLEREGISTRY_H
#define GEL_VARIABLEREGISTRY_H

#include "Console.h"

namespace GLGraphics {

    namespace Registry {
        class command {
        public:
            inline command() : m_console(NULL) {}

            inline void reg(Console &cs,
                            const std::string &name,
                            const std::function<void()> &function,
                            const std::string &help) {
                assert(!m_console);
                m_console = &cs;
                m_id = m_console->reg_cmd0(name, function, help);
            }

            template<typename A0>
            void reg(Console &cs,
                     const std::string &name,
                     const std::function<void(const A0 &)> &function,
                     const std::string &help) {
                assert(!m_console);
                m_console = &cs;
                m_id = m_console->reg_cmd1<A0>(name, function, help);
            }

            template<typename A0, typename A1>
            void reg(Console &cs,
                     const std::string &name,
                     const std::function<void(const A0 &, const A1 &)> &function,
                     const std::string &help) {
                assert(!m_console);
                m_console = &cs;
                m_id = m_console->reg_cmd2<A0, A1>(name, function, help);
            }

            template<typename A0, typename A1, typename A2>
            void reg(Console &cs,
                     const std::string &name,
                     const std::function<void(const A0 &,
                                              const A1 &, const A2 &)> &function,
                     const std::string &help) {
                assert(!m_console);
                m_console = &cs;
                m_id = m_console->reg_cmd3<A0, A1, A2>(name, function, help);
            }

            inline ~command() {
                if (m_console)
                    m_console->unreg_cmd(m_id);
            }

            inline const char *get_name() const {
                assert(m_console);
                return m_console->get_name(m_id);
            }

            inline const char *get_help() const {
                assert(m_console);
                return m_console->get_help(m_id);
            }

            inline Console *get_console() const { return m_console; }

            inline Console::cmd_token get_id() const {
                assert(m_console);
                return m_id;
            }

        private:
            command(command &);

            const command &operator=(const command &);

            Console *m_console;
            Console::cmd_token m_id;
        };

        template<typename T>
        class variable {
        public:
            variable(const T &initial_value = T())
                    : m_value(initial_value) {}

            void reg(Console &cs,
                     const std::string &name,
                     const std::string &help) {
                if (m_set_cmd.get_console() == 0) {
                    m_print_cmd.reg(cs, name,
                                    std::bind(&variable::print_value, this), help);

                    m_set_cmd.reg<T>(cs, name,
                                     std::bind(&variable::set_value, this, std::placeholders::_1),
                                     help);
                }
            }

            const variable &operator=(const T &value) {
                m_value = value;
                return *this;
            }

            operator const T &() const { return m_value; }

            const char *get_name() const { return m_print_cmd.get_name(); }

            const char *get_help() const { return m_print_cmd.get_help(); }

        private:
            variable(const variable &);

            const variable &operator=(const variable &);

            void print_value() {
                m_print_cmd.get_console()->printf("%s = %s",
                                                  m_print_cmd.get_name(),
                                                  Console::to_string(m_value).c_str());
            }

            void set_value(const T &value) {
                m_value = value;
                m_print_cmd.get_console()->execute(m_print_cmd.get_name());
            }

            T m_value;

            command m_print_cmd;
            command m_set_cmd;
        };
    }
}
#endif //GEL_VARIABLEREGISTRY_H
