import { Position, Toaster, Intent } from "@blueprintjs/core";

/** Singleton toaster instance. Create separate instances for different options. */

const ToastTopCenter = Toaster.create({
  className: "recipe-toaster",
  position: Position.TOP
});

const ToastBottomCenter = Toaster.create({
  className: "recipe-toaster",
  position: Position.BOTTOM
});

/*
A "user" error - eg, bad input
*/
export const postUserErrorToast = message =>
  ToastTopCenter.show({ message, intent: Intent.WARNING });

/*
A toast the user must dismiss manually, because they need to act on its information,
ie., 8 bulk add genes out of 40 were bad. Manually see which ones and fix.
*/
export const keepAroundErrorToast = message =>
  ToastTopCenter.show({ message, timeout: 0, intent: Intent.WARNING });

/*
a hard network error
*/
export const postNetworkErrorToast = message =>
  ToastTopCenter.show({
    message,
    timeout: 30000,
    intent: Intent.DANGER
  });

/*
Async message to user
*/
export const postAsyncSuccessToast = message =>
  ToastTopCenter.show({
    message,
    timeout: 10000,
    intent: Intent.SUCCESS
  });

export const postAsyncFailureToast = message =>
  ToastTopCenter.show({
    message,
    timeout: 10000,
    intent: Intent.WARNING
  });

export const termsOfServiceToast = (message, _onDismiss) =>
  ToastBottomCenter.show({
    message,
    timeout: 0,
    intent: Intent.PRIMARY,
    onDismiss: _onDismiss
  });
