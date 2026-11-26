/*
 *    ViennaRNA/utils/basic.c
 *
 *               c  Ivo L Hofacker and Walter Fontana
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* isatty() is not available for Windows systems */
#ifndef _WIN32
#include <unistd.h>
#endif

#if VRNA_WITH_PTHREADS
# include <pthread.h>
#endif

#include "ViennaRNA/intern/color_output.h"

#include "ViennaRNA/utils/basic.h"

#define EXIT_ON_ERROR

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
typedef struct {
  vrna_log_cb_f       cb;
  void                *cb_data;
  vrna_logdata_free_f data_release;
  vrna_log_levels_e   level;
} logger_callback;

PRIVATE struct {
  FILE              *default_file;
  vrna_log_levels_e default_level;
  unsigned int      options;
  vrna_log_lock_f   lock;
  void              *lock_data;
  vrna_array(logger_callback) callbacks;
#if VRNA_WITH_PTHREADS
  pthread_mutex_t   mtx;            /* semaphore to prevent concurrent access */
#endif
} logger = {
  .default_file   = NULL,
  .default_level  = VRNA_LOG_LEVEL_DEFAULT,
  .options        = VRNA_LOG_OPTION_DEFAULT,
  .lock           = NULL,
  .lock_data      = NULL,
  .callbacks      = NULL,
#if VRNA_WITH_PTHREADS
  .mtx            = PTHREAD_MUTEX_INITIALIZER
#endif
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
log_v(vrna_log_event_t *event);

PRIVATE void
log_default(vrna_log_event_t *event);

PRIVATE void
lock(void);

PRIVATE void
unlock(void);

#ifndef VRNA_WITHOUT_TTY_COLORS
PRIVATE const char *
get_log_level_color(vrna_log_levels_e level);
#endif

PRIVATE const char *
get_log_level_string(vrna_log_levels_e level);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...)
{
  vrna_log_event_t event = {
    .format_string  = format_string,
    .level          = level,
    .line_number    = line_number,
    .file_name      = file_name
  };

  va_start(event.params, format_string);
  log_v(&event);
  va_end(event.params);
}

/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
log_default(vrna_log_event_t *event)
{
  if (!logger.default_file)
    logger.default_file = stderr;

  /* print time unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_TIME) {
    char    timebuf[64];
    time_t  t = time(NULL);
    timebuf[strftime(timebuf, sizeof(timebuf), "%H:%M:%S", localtime(&t))] = '\0';
    fprintf(logger.default_file, "%s ", timebuf);
  }

  /* print log level */
#ifndef VRNA_WITHOUT_TTY_COLORS
  if (isatty(fileno(logger.default_file))) {
    fprintf(logger.default_file,
            "%s%-9s" ANSI_COLOR_RESET " ",
            get_log_level_color(event->level),
            get_log_level_string(event->level));
  } else {
#endif
  fprintf(logger.default_file,
          "%-9s ",
          get_log_level_string(event->level));
#ifndef VRNA_WITHOUT_TTY_COLORS
}

#endif

  /* print file name / line number trace unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_CALL) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (isatty(fileno(logger.default_file))) {
      fprintf(logger.default_file,
              "\x1b[90m%s:%d:" ANSI_COLOR_RESET " ",
              event->file_name,
              event->line_number);
    } else {
#endif
    fprintf(logger.default_file, "%s:%d: ",
            event->file_name,
            event->line_number);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  /* print actual message */
  vfprintf(logger.default_file, event->format_string, event->params);
  fprintf(logger.default_file, "\n");
  fflush(logger.default_file);
}

PRIVATE void
log_v(vrna_log_event_t *event)
{
  lock();

  /* initialize the logger, if not done already */
  if (logger.callbacks == NULL) {
    /* initialize callbacks if not done so far */
    vrna_array_init(logger.callbacks);
  }

  /* process log output for default implementation */
  if (!(logger.options & VRNA_LOG_OPTION_QUIET)) {
    /* print log if not in quiet mode */
    if (event->level >= logger.default_level)
      log_default(event);
  }

  /* process log for any user-defined output */
  for (size_t i = 0; i < vrna_array_size(logger.callbacks); i++) {
    logger_callback *cb = &(logger.callbacks[i]);

    if (event->level >= cb->level)
      cb->cb(event, cb->cb_data);
  }

  unlock();
}

PRIVATE void
lock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(1, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_lock(&(logger.mtx));
}
#endif
}

PRIVATE void
unlock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(0, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_unlock(&(logger.mtx));
}
#endif
}

PRIVATE const char *
get_log_level_string(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return "[DEBUG]";
    case VRNA_LOG_LEVEL_INFO:
      return "[INFO]";
    case VRNA_LOG_LEVEL_WARNING:
      return "[WARNING]";
    case VRNA_LOG_LEVEL_ERROR:
      return "[ERROR]";
    case VRNA_LOG_LEVEL_CRITICAL:
      return "[FATAL]";
    default:
      return "[UNKNOWN]";
  }
}

#ifndef VRNA_WITHOUT_TTY_COLORS
PRIVATE const char *
get_log_level_color(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return ANSI_COLOR_CYAN_B;
    case VRNA_LOG_LEVEL_INFO:
      return ANSI_COLOR_BLUE_B;
    case VRNA_LOG_LEVEL_WARNING:
      return ANSI_COLOR_YELLOW_B;
    case VRNA_LOG_LEVEL_ERROR:
      return ANSI_COLOR_RED_B;
    case VRNA_LOG_LEVEL_CRITICAL:
      return ANSI_COLOR_MAGENTA_B;
    default:
      return ANSI_COLOR_RESET;
  }
}
#endif

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#endif
