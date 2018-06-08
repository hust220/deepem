#pragma once

#include <functional>
#include <vector>
#include <string>
#include <iostream>
#include <csignal>
#include <csetjmp>

struct BaseCase
{
    virtual void run() = 0;
    virtual ~BaseCase() = default;
};

class UnitTest
{
    private:
        UnitTest() : last_checked_line_{ 0 }, failure_num_{ 0 } {}
        ::std::vector<BaseCase*> test_cases_;
        ::std::string last_checked_file_;
        size_t last_checked_line_;
        size_t failure_num_;

    public:
        static jmp_buf& abort_case()
        {
            static jmp_buf buf;
            return buf;
        }

        static UnitTest& getInstance()
        {
            static UnitTest instance;
            return instance;
        }

        void runAll()
        {
            ::std::cout << ">>> running " << test_cases_.size() << " tests..." << ::std::endl;
            for (BaseCase *test : test_cases_)
            {
                if (setjmp(abort_case()))
                    continue;
                test->run();
            }
        }

        void registerTestCase(BaseCase *test)
        {
            test_cases_.push_back(test);
        }

        void printLastCheckedPoint()
        {
            ::std::cout << ">>> ";
            ::std::cout << last_checked_file_ << "(" << last_checked_line_ << ")" << ": last checkpoint" << ::std::endl;
        }

        void checkFile(const ::std::string& file)
        {
            last_checked_file_ = file;
        }

        void checkLine(size_t line)
        {
            last_checked_line_ = line;
        }

        void incFailure()
        {
            ++failure_num_;
        }

        size_t getFailureNum()
        {
            return failure_num_;
        }
};

struct TestCase : BaseCase
{
    public:
        TestCase(::std::function<void()> method, const ::std::string& name, const ::std::string& file, size_t line)
            : method_{ method }, case_name_{ name }, defined_file_{ file }, defined_line_{ line }
        {
            UnitTest::getInstance().registerTestCase(this);
        }

        void run() override
        {
#ifdef NDEBUG
            try
            {
#endif
                UnitTest::getInstance().checkFile(defined_file_);
                UnitTest::getInstance().checkLine(defined_line_);
                size_t old_failure_num = UnitTest::getInstance().getFailureNum();
                method_();
                int failures = UnitTest::getInstance().getFailureNum() - old_failure_num;
                if (failures)
                {
                    ::std::cout << ">>> ";
                    ::std::cout << failures << " failures are detected in the test case \"" << case_name_ << "\"" << ::std::endl;
                }
#ifdef NDEBUG
            }
            catch (::std::exception& e)
            {
                UnitTest::getInstance().incFailure();
                ::std::cout << ">>> fatal error: in \"" << case_name_ << "\": " << typeid(e).name() << ": " << e.what() << ::std::endl;
                UnitTest::getInstance().printLastCheckedPoint();
            }
            catch (...)
            {
                UnitTest::getInstance().incFailure();
                ::std::cout << ">>> fatal error: in \"" << case_name_ << "\": unknown type exception" << ::std::endl;
                UnitTest::getInstance().printLastCheckedPoint();
            }
#endif
        }

        ~TestCase() override = default;

    private:
        ::std::function<void()> method_;
        ::std::string case_name_;
        ::std::string defined_file_;
        size_t defined_line_;
};

    [[noreturn]]
static void report_and_exit()
{
    ::std::cout << "\n**** ";
    ::std::cout << UnitTest::getInstance().getFailureNum() << " failures are detected." << ::std::endl;
    exit(UnitTest::getInstance().getFailureNum());
}

namespace unit_test {

    static void run_test()
    {
#ifdef NDEBUG
        signal(SIGSEGV, [](int)
                {
                UnitTest::getInstance().incFailure();
                ::std::cout << ">>> fatal error: received SIGSEGV." << ::std::endl;
                UnitTest::getInstance().printLastCheckedPoint();
                report_and_exit();
                });
#endif
        UnitTest::getInstance().runAll();
        report_and_exit();
    }

}

    template <typename F, typename... Args, typename = decltype(::std::declval<F>()(::std::declval<Args>()...))>
void do_check_failed(F&& f, Args&&... args)
{
    f(::std::forward<Args>(args)...);
}

    template <typename... Msgs>
void do_check_failed(Msgs&&... msgs)
{
    (void)::std::initializer_list<int>{(::std::cout << ">>> " << msgs << ::std::endl, 0)...};
}

#define TEST_CASE(test_name)                                                                        \
    static void test_name();                                                                        \
    static TestCase test_name##_case{test_name, #test_name, __FILE__, __LINE__};                    \
    static void test_name()

#define G_CHECK(cond, strict, ...)                                                                  \
    do {                                                                                            \
        UnitTest::getInstance().checkFile(__FILE__);                                                \
        UnitTest::getInstance().checkLine(__LINE__);                                                \
        if(!(cond)) {                                                                               \
            UnitTest::getInstance().incFailure();                                                   \
            if(strict) {                                                                            \
                ::std::cout << ">>> check \"" << #cond << "\" failed." << ::std::endl;                  \
                ::std::cout << ">>> critical error at " __FILE__ "(" << __LINE__ << ")." << ::std::endl;\
                do_check_failed(__VA_ARGS__);                                                       \
                longjmp(UnitTest::abort_case(), true);                                              \
            } else {                                                                                \
                ::std::cout << ">>> check \"" << #cond << "\" failed." << "at "                       \
                << __FILE__ << "(" << __LINE__ << ")" << ::std::endl;                             \
                do_check_failed(__VA_ARGS__);                                                       \
            }                                                                                       \
        }                                                                                           \
    } while(0)

#define TEST_CHECK(cond, ...)                                                                   \
    G_CHECK(cond, false, __VA_ARGS__)

#define TEST_REQUIRE(cond, ...)                                                                 \
    G_CHECK(cond, true, __VA_ARGS__)
